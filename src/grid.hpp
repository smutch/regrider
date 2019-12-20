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

#ifndef GRID_H
#define GRID_H

#include <array>
#include <complex>
#include <memory>

/** A 3D grid class to handle input independent functionality.
 */
class Grid {
public:
    std::array<int32_t, 3> n_cell; /** The number of cells in each dimension */
    std::array<double, 3> box_size; /** The box size in input units (typically h^-1 Mpc) */
    int n_logical; /** The total number of cells */
    int n_padded; /** Number of elements in the padded array */
    int n_complex; /** The number of complex elements in the FFTd array */
    bool flag_padded = false; /** Has the indexing been reorder to be padded for an inplace FFT? */

private:
    std::unique_ptr<float, void (*)(float*)> grid; /** A pointer to the grid data, allowing it to be
                                                      automatically freed when this Grid object goes out
                                                      of scope. */

public:
    /** The indexing type required for an ::index function call.
     */
    enum class index_type {
        padded,
        real,
        complex_herm
    };

    /** The type of filter to be applied to the grid to smooth it.
     */
    enum class filter_type {
        real_top_hat,
        k_top_hat,
        gaussian
    };

    /** Basic constructor.
     * This will allocate the grid array, and store the corresponding size in various forms.
     *
     * @param n_cell_ The number of @logical cells in each dimension
     * @param box_size_ The size of the simulation volume in input units
     */
    Grid(const std::array<int32_t, 3> n_cell_, const std::array<double, 3> box_size_);

    /** Update the "size" of the grid for a new logical size.
     * Note that this does not alter the size of the memory allocation, just what this allocation represents.
     *
     * @param n_cell_ The new number of @logical cells in each dimension
     */
    void update_properties(const std::array<int32_t, 3> n_cell_);

    /** Return the pointer to the grid data.
     *
     * @return Float pointer to the grid data.
     */
    float* get();

    /** Return the pointer to the grid data, cast as a complex array.
     *
     * @return Complex pointer to the grid data
     */
    std::complex<float>* get_complex();

    /** Indexing function for arbitrary grid of any 3D size.
     *
     * @param i Index in 1st dimension
     * @param j Index in second dimension
     * @param k Index in third dimension
     * @param shape Shape of the 3D array
     * @return The index
     */
    constexpr int index(int i, int j, int k, index_type type, std::array<int, 3> shape);

    /** Indexing function for the current grid.
     *
     * @param i Index in 1st dimension
     * @param j Index in second dimension
     * @param k Index in third dimension
     * @return The index
     */
    constexpr int index(int i, int j, int k, index_type type);

    /** Convert the grid from logical memory ordering to padded ordering.
     */
    void real_to_padded_order();

    /** Convert the grid from padded memory ordering to logical ordering.
     */
    void padded_to_real_order();

    /** Do the forward FFT
     *
     * @param n_threads The number of openmp threads ot use for the transform. -1 --> all available (default).
     */
    void forward_fft(int n_threads = -1);

    /** Do the reverse FFT
     *
     * @param n_threads The number of openmp threads ot use for the transform. -1 --> all available (default).
     */
    void reverse_fft(int n_threads = -1);

    /** Filter the grid using a given filter type and size.
     *
     * @param type The filter type to use
     * @param R the size (typically radius) of the filter
     */
    void filter(filter_type type, const float R);

    /** Subsample the grid to provide a new one with the requested dimensions.
     *
     * Note that the parameters of the Grid object will be updated correspondingly.
     *
     * @param new_n_cell The new logical size of the grid.
     */
    void sample(const std::array<int, 3> new_n_cell);
};

#endif
