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
#include <fftw3.h>
#include <omp.h>


class Grid {
    public:
    std::array<int32_t, 3> n_cell;
    int n_logical;
    int n_padded;
    int n_complex;

    Grid(std::array<int32_t, 3>n_cell_) :
        n_cell{n_cell_},
        n_logical{n_cell[0] * n_cell[1] * n_cell[2]},
        n_padded{n_cell[0] * n_cell[1] * 2*(n_cell[2]/2+1)}, 
        n_complex{n_cell[0] * n_cell[1] * (n_cell[2]/2+1)},
        grid(fftwf_alloc_real(n_padded), [](float* grid){ fftwf_free(grid); })
        {};

    float* get() {
        return grid.get();
    }

    fftwf_complex* get_complex() {
        return (fftwf_complex*)grid.get();
    }

    enum index_type {
        padded,
        real
    };

    constexpr int index(int i, int j, int k, index_type type) {
        auto index = k + n_cell[1] * (j + n_cell[0] * i);
        assert(index < n_logical);

        switch (type) {
            case padded:
                index = k + (2 * (n_cell[1] / 2 + 1)) * (j + n_cell[0] * i);
                break;
            case real:
                break;
            default:
                fmt::print(stderr, "Unrecognised index_type!\n");
                break;
        }

        return index;
    }

    void real_to_padded_order() {
        for (int32_t ii = n_cell[0]-1; ii >= 0; --ii)
            for (int32_t jj = n_cell[1] - 1; jj >= 0; --jj)
                for (int32_t kk = n_cell[2] - 1; kk >= 0; --kk)
                    grid.get()[index(ii, jj, kk, index_type::padded)] = grid.get()[index(ii, jj, kk, index_type::real)];
    }

    void padded_to_real_order() {
        for (int32_t ii = n_cell[0]-1; ii >= 0; --ii)
            for (int32_t jj = n_cell[1] - 1; jj >= 0; --jj)
                for (int32_t kk = n_cell[2] - 1; kk >= 0; --kk)
                    grid.get()[index(ii, jj, kk, index_type::real)] = grid.get()[index(ii, jj, kk, index_type::padded)];
    }

    private:
    std::unique_ptr<float, void(*)(float*)> grid;
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

    auto orig = Grid(n_cell);

    bool found = false;
    for(int32_t ii=0; ii<n_grids && !found; ++ii){

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

    orig.real_to_padded_order();

    fftwf_plan_with_nthreads(omp_get_max_threads());
    auto plan = fftwf_plan_dft_r2c_3d(n_cell[0], n_cell[1], n_cell[2], orig.get(), orig.get_complex(), FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    // Remember to add the factor of VOLUME/TOT_NUM_PIXELS when converting from
    // real space to k-space.
    // Note: we will leave off factor of VOLUME, in anticipation of the inverse
    // FFT below
    for (int ii = 0; ii < orig.n_complex; ++ii)
        orig.get_complex()[ii][0] /= orig.n_logical;

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
