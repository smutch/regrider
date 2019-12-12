/*
 * regrider: Downsample gbpTrees and VELOCIraptor grids using FFTW
 * Copyright © 2019 Simon Mutch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <cxxopts.hpp>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cstring>
#include <complex>
#include <fftw3.h>
#include <omp.h>


class Grid {
    public:
    std::array<int32_t, 3> n_cell;
    std::array<double, 3> box_size;
    int n_logical;
    int n_padded;
    int n_complex;
    bool flag_padded = false;

    private:
    std::unique_ptr<float, void(*)(float*)> grid;

    public:
    enum index_type {
        padded,
        real,
        complex_herm
    };
    
    enum filter_type {
        real_top_hat,
        k_top_hat,
        gaussian
    };

    Grid(const std::array<int32_t, 3>n_cell_, const std::array<double, 3>box_size_) :
        n_cell{n_cell_},
        box_size{box_size_},
        n_logical{n_cell[0] * n_cell[1] * n_cell[2]},
        n_padded{n_cell[0] * n_cell[1] * 2*(n_cell[2]/2+1)}, 
        n_complex{n_cell[0] * n_cell[1] * (n_cell[2]/2+1)},
        grid(fftwf_alloc_real(n_padded), [](float* grid){ fftwf_free(grid); })
        {};

    float* get() {
        return grid.get();
    }

    std::complex<float>* get_complex() {
        return (std::complex<float>*)grid.get();
    }

    constexpr int index(int i, int j, int k, index_type type) {
        auto index = k + n_cell[1] * (j + n_cell[0] * i);
        assert(index < n_logical);

        switch (type) {
            case padded:
                index = k + (2 * (n_cell[1] / 2 + 1)) * (j + n_cell[0] * i);
                break;
            case real:
                break;
            case complex_herm:
                index = k + (n_cell[1] / 2 + 1) * (j + n_cell[0] * i);
                break;
            default:
                fmt::print(stderr, "Unrecognised index_type!\n");
                break;
        }

        return index;
    }

    void real_to_padded_order() {
        for (int ii = n_cell[0]-1; ii >= 0; --ii)
            for (int jj = n_cell[1] - 1; jj >= 0; --jj)
                for (int kk = n_cell[2] - 1; kk >= 0; --kk)
                    grid.get()[index(ii, jj, kk, index_type::padded)] = grid.get()[index(ii, jj, kk, index_type::real)];
        flag_padded = true;
    }

    void padded_to_real_order() {
        for (int ii = n_cell[0]-1; ii >= 0; --ii)
            for (int jj = n_cell[1] - 1; jj >= 0; --jj)
                for (int kk = n_cell[2] - 1; kk >= 0; --kk)
                    grid.get()[index(ii, jj, kk, index_type::real)] = grid.get()[index(ii, jj, kk, index_type::padded)];
        flag_padded = false;
    }

    void forward_fft(int n_threads = -1) {
        if (n_threads == -1) {
            n_threads = omp_get_max_threads();
        }

        if (!flag_padded) {
            real_to_padded_order();
        }

        fftwf_plan_with_nthreads(n_threads);
        auto plan = fftwf_plan_dft_r2c_3d(n_cell[0], n_cell[1], n_cell[2], get(), (fftwf_complex*)get_complex(), FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        // Remember to multiply by VOLUME/TOT_NUM_PIXELS when converting from
        // real space to k-space.  Note: we will leave off factor of VOLUME, in
        // anticipation of the inverse FFT
        for (int ii = 0; ii < n_complex; ++ii)
            get_complex()[ii] /= n_logical;
    }

    void reverse_fft(int n_threads = -1) {
        if (n_threads == -1) {
            n_threads = omp_get_max_threads();
        }

        fftwf_plan_with_nthreads(n_threads);
        auto plan = fftwf_plan_dft_c2r_3d(n_cell[0], n_cell[1], n_cell[2], (fftwf_complex*)get_complex(), get(), FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        padded_to_real_order();
    }

    void filter(filter_type type, const float R) {

        forward_fft();

        const int middle = n_cell[2] / 2;
        std::array<double, 3> delta_k = {0};

        for (int ii=0; ii < 3; ++ii) {
            delta_k[ii] = (2.0 * M_PI / box_size[ii]);
        }

        // Loop through k-box
        for (int n_x = 0; n_x < n_cell[0]; ++n_x) {
            double k_x = 0;

            if (n_x > middle)
                k_x = (n_x - n_cell[0]) * delta_k[0];
            else
                k_x = n_x * delta_k[0];

            for (int n_y = 0; n_y < n_cell[1]; ++n_y) {
                double k_y = 0;

                if (n_y > middle)
                    k_y = (n_y - n_cell[1]) * delta_k[1];
                else
                    k_y = n_y * delta_k[1];

                for (int n_z = 0; n_z <= middle; ++n_z) {
                    double k_z = n_z * delta_k[2];

                    double k_mag = sqrt(k_x * k_x + k_y * k_y + k_z * k_z);

                    double kR = k_mag * R;

                    switch (type) {
                        case real_top_hat: // Real space top-hat
                            if (kR > 1e-4) {
                                get_complex()[index(n_x, n_y, n_z, complex_herm)] *= 3.0 * (sin(kR) / pow(kR, 3) - cos(kR) / pow(kR, 2));
                            }
                            break;

                        case k_top_hat: // k-space top hat
                            kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
                            if (kR > 1) {
                                get_complex()[index(n_x, n_y, n_z, complex_herm)] = 0.0;
                            }
                            break;

                        case gaussian: // Gaussian
                            kR *= 0.643; // Equates integrated volume to the real space top-hat
                            get_complex()[index(n_x, n_y, n_z, complex_herm)] *= pow(M_E, -kR * kR / 2.0);
                            break;

                        default:
                            if ((n_x == 0) && (n_y == 0) && (n_z == 0)) {
                                fmt::print(stderr, "Error: filter type {} is undefined!", type);
                            }
                            break;
                    }
                }
            }
        } // End looping through k box

        reverse_fft();
    }

    void sample(const std::array<int, 3> new_n_cell) {

        std::array<int, 3> n_every = {0};
        for(int ii=0; ii < 3; ++ii) {
            n_every[ii] = n_cell[ii] / new_n_cell[ii];
        }

        // TODO cont here...

    }

};


static void read_gbptrees(const std::string fname_in, const std::string grid_name)
{
    fmt::print("Reading gbpTrees grid file {}...\n", fname_in);
    std::ifstream ifs(fname_in, std::ios::binary | std::ios::in);

    std::array<int, 3> n_cell;
    ifs.read((char*)(n_cell.data()), sizeof(int) * 3);
    fmt::print("n_cell = {}\n", fmt::join(n_cell, ","));

    std::array<double, 3> box_size;
    ifs.read((char*)(box_size.data()), sizeof(double) * 3);
    fmt::print("box_size = {}\n", fmt::join(box_size, ","));

    int32_t n_grids;
    ifs.read((char*)(&n_grids), sizeof(int));
    fmt::print("n_grids = {}\n", n_grids);

    int32_t ma_scheme;
    ifs.read((char*)(&ma_scheme), sizeof(int));
    fmt::print("ma_scheme = {}\n", ma_scheme);

    auto orig = Grid(n_cell, box_size);

    bool found = false;
    for(int ii=0; ii<n_grids && !found; ++ii){

        std::string ident(32, '\0');
        ifs.read((char*)(ident.data()), sizeof(ident));
        ident.resize(strlen(ident.c_str()));
        
        if (grid_name == ident) {
            fmt::print("Reading grid {}...\n", ident);
            ifs.read((char*)orig.get(), sizeof(float)*orig.n_logical);
            found = true;
        } else {
            fmt::print("Skipping grid {}...\n", ident);
            ifs.seekg(sizeof(float)*orig.n_logical, std::ifstream::cur);
        }

    }

    ifs.close();

    if (!found) {
        fmt::print(stderr, "Failed to find grid named `{}' in file!\n", grid_name);
        return;
    }

    // DEBUG
    {
        std::vector subset(orig.get(), orig.get() + 10);
        fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }

    orig.forward_fft();

    // filter(slab,
            // (int)slab_ix_start,
            // (int)slab_nix, n_cell[0],
            // (float)(run_globals.params.BoxSize / (double)run_globals.params.ReionGridDim / 2.0));

    fmt::print("...done\n");
}

static void write_gbptrees(const std::string fname_out) {}

int main(int argc, char* argv[])
{
    cxxopts::Options options(
        "regrider", "Downsample gbpTrees and VELOCIraptor trees using FFTW");

    options.add_options()("d,dim", "new grid dimension", cxxopts::value<uint32_t>())
        ("n,name", "grid name (must match conventions), ", cxxopts::value<std::string>())
        ("g,gbptrees", "input gbpTrees grid file", cxxopts::value<std::string>())
        ("v,velociraptor", "input VELOCIraptor grid file", cxxopts::value<std::string>());

    auto vm = options.parse(argc, argv);

    if (vm.count("gbptrees") && vm.count("velociraptor")) {
        fmt::print(stderr, "Must specify either gbpTrees or VELOCIraptor file. Not both...\n");
        return 1;
    }

    fftwf_init_threads();

    if (vm.count("gbptrees")) {
        read_gbptrees(vm["gbptrees"].as<std::string>(), vm["name"].as<std::string>());
    }

    fftwf_cleanup_threads();

    return 0;
}
