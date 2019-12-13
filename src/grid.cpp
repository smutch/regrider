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
#include <fftw3.h>
#include <cassert>
#include <omp.h>
#include <iostream>
#include <vector>

#include "grid.hpp"

Grid::Grid(const std::array<int32_t, 3> n_cell_, const std::array<double, 3> box_size_)
    : n_cell { n_cell_ }
    , box_size { box_size_ }
    , n_logical { n_cell[0] * n_cell[1] * n_cell[2] }
    , n_padded { n_cell[0] * n_cell[1] * 2 * (n_cell[2] / 2 + 1) }
    , n_complex { n_cell[0] * n_cell[1] * (n_cell[2] / 2 + 1) }
    , grid(fftwf_alloc_real(n_padded), [](float* grid) { fftwf_free(grid); }) {};

void Grid::update_properties(const std::array<int32_t, 3> n_cell_) {
    n_cell = n_cell_;
    n_logical = n_cell[0] * n_cell[1] * n_cell[2];
    n_padded = n_cell[0] * n_cell[1] * 2 * (n_cell[2] / 2 + 1);
    n_complex = n_cell[0] * n_cell[1] * (n_cell[2] / 2 + 1);
}

float* Grid::get()
{
    return grid.get();
}

std::complex<float>* Grid::get_complex()
{
    return (std::complex<float>*)grid.get();
}

constexpr int Grid::index(int i, int j, int k, index_type type, std::array<int, 3>shape)
{
    auto index = k + shape[1] * (j + shape[0] * i);
    assert(index < n_logical);

    switch (type) {
    case index_type::padded:
        index = k + (2 * (shape[1] / 2 + 1)) * (j + shape[0] * i);
        break;
    case index_type::real:
        break;
    case index_type::complex_herm:
        index = k + (shape[1] / 2 + 1) * (j + shape[0] * i);
        break;
    default:
        fmt::print(stderr, "Unrecognised index_type!\n");
        break;
    }

    return index;
}

constexpr int Grid::index(int i, int j, int k, index_type type)
{
    return Grid::index(i, j, k, type, n_cell);
}

void Grid::real_to_padded_order()
{
    auto grid_ = grid.get();
    for (int ii = n_cell[0] - 1; ii >= 0; --ii)
        for (int jj = n_cell[1] - 1; jj >= 0; --jj)
            for (int kk = n_cell[2] - 1; kk >= 0; --kk)
                grid_[index(ii, jj, kk, index_type::padded)] = grid_[index(ii, jj, kk, index_type::real)];
    flag_padded = true;
}

void Grid::padded_to_real_order()
{
    auto grid_ = get();
    for (int ii = n_cell[0] - 1; ii >= 0; --ii)
        for (int jj = n_cell[1] - 1; jj >= 0; --jj)
            for (int kk = n_cell[2] - 1; kk >= 0; --kk)
                grid_[index(ii, jj, kk, index_type::real)] = grid_[index(ii, jj, kk, index_type::padded)];
    flag_padded = false;
}

void Grid::forward_fft(int n_threads)
{
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
    auto complex_grid = get();
#pragma omp parallel for default(none) private(n_logical) shared(complex_grid)
    for (int ii = 0; ii < n_complex; ++ii)
        complex_grid[ii] /= n_logical;
}

void Grid::reverse_fft(int n_threads)
{
    if (n_threads == -1) {
        n_threads = omp_get_max_threads();
    }

    fftwf_plan_with_nthreads(n_threads);
    auto plan = fftwf_plan_dft_c2r_3d(n_cell[0], n_cell[1], n_cell[2], (fftwf_complex*)get_complex(), get(), FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    padded_to_real_order();
}

void Grid::filter(filter_type type, const float R)
{

    fmt::print("Filtering grid: ");
    std::cout << std::flush;

    fmt::print("doing forward fft... ");
    std::cout << std::flush;

    forward_fft();

    const int middle = n_cell[2] / 2;
    std::array<double, 3> delta_k = { 0 };

    for (int ii = 0; ii < 3; ++ii) {
        delta_k[ii] = (2.0 * M_PI / box_size[ii]);
    }

    // Loop through k-box
    fmt::print("applying convolution... ");
    std::cout << std::flush;

    auto complex_grid = get_complex();

#pragma omp parallel for default(none) private(delta_k) shared(complex_grid, stderr, type)
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
                case filter_type::real_top_hat: // Real space top-hat
                    if (kR > 1e-4) {
                        complex_grid[index(n_x, n_y, n_z, index_type::complex_herm)] *= 3.0 * (sin(kR) / pow(kR, 3) - cos(kR) / pow(kR, 2));
                    }
                    break;

                case filter_type::k_top_hat: // k-space top hat
                    kR *= 0.413566994; // Equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
                    if (kR > 1) {
                        complex_grid[index(n_x, n_y, n_z, index_type::complex_herm)] = 0.0;
                    }
                    break;

                case filter_type::gaussian: // Gaussian
                    kR *= 0.643; // Equates integrated volume to the real space top-hat
                    complex_grid[index(n_x, n_y, n_z, index_type::complex_herm)] *= pow(M_E, -kR * kR / 2.0);
                    break;

                default:
                    if ((n_x == 0) && (n_y == 0) && (n_z == 0)) {
                        fmt::print(stderr, "Error: filter type {} is undefined!", (int)type);
                    }
                    break;
                }
            }
        }
    } // End looping through k box

    fmt::print("doing inverse fft... ");
    std::cout << std::flush;

    reverse_fft();

    fmt::print("done\n");
}

void Grid::sample(const std::array<int, 3> new_n_cell)
{
    fmt::print("Subsampling grid... ");

    std::array<int, 3> n_every = { 0 };
    for (int ii = 0; ii < 3; ++ii) {
        n_every[ii] = n_cell[ii] / new_n_cell[ii];
    }

    auto grid_ = grid.get();

    // TODO: I need to check to make sure this is valid
    for(int ii = 0, ii_lo = 0; ii < n_cell[0]; ii += n_every[0], ++ii_lo) {
        for(int jj = 0, jj_lo = 0; jj < n_cell[1]; jj += n_every[1], ++jj_lo) {
            for(int kk = 0, kk_lo = 0; kk < n_cell[2]; kk += n_every[2], ++kk_lo) {
                grid_[index(ii_lo, jj_lo, kk_lo, index_type::real, new_n_cell)] = grid_[index(ii, jj, kk, index_type::real)];
            }
        }
    }

    update_properties(new_n_cell);

    fmt::print("done\n");
}
