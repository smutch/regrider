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

#include <array>
#include <memory>
#include <complex>

class Grid {
public:
    std::array<int32_t, 3> n_cell;
    std::array<double, 3> box_size;
    int n_logical;
    int n_padded;
    int n_complex;
    bool flag_padded = false;

private:
    std::unique_ptr<float, void (*)(float*)> grid;

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

    Grid(const std::array<int32_t, 3> n_cell_, const std::array<double, 3> box_size_);

    float* get();
    std::complex<float>* get_complex();
    constexpr int index(int i, int j, int k, index_type type, std::array<int, 3>shape);
    constexpr int index(int i, int j, int k, index_type type);
    void real_to_padded_order();
    void padded_to_real_order();
    void forward_fft(int n_threads = -1);
    void reverse_fft(int n_threads = -1);
    void filter(filter_type type, const float R);
    void sample(const std::array<int, 3> new_n_cell);
};
