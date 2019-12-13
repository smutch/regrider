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
#include <vector>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fstream>

#include "gbptrees.hpp"

void read_gbptrees(const std::string fname_in, const std::string grid_name, Grid& grid)
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

    grid.init(n_cell, box_size);

    bool found = false;
    for (int ii = 0; ii < n_grids && !found; ++ii) {

        std::string ident(32, '\0');
        ifs.read((char*)(ident.data()), sizeof(ident));
        ident.resize(strlen(ident.c_str()));

        if (grid_name == ident) {
            fmt::print("Reading grid {}... ", ident);
            ifs.read((char*)grid.get(), sizeof(float) * grid.n_logical);
            found = true;
        } else {
            fmt::print("Skipping grid {}... ", ident);
            ifs.seekg(sizeof(float) * grid.n_logical, std::ifstream::cur);
        }
    }

    ifs.close();

    fmt::print("done\n");

    if (!found) {
        fmt::print(stderr, "Failed to find grid named `{}' in file!\n", grid_name);
        exit(EXIT_FAILURE);
    }

    // DEBUG
    {
        std::vector subset(grid.get(), grid.get() + 10);
        fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }

    // filter(slab,
    // (int)slab_ix_start,
    // (int)slab_nix, n_cell[0],
    // (float)(run_globals.params.BoxSize / (double)run_globals.params.ReionGridDim / 2.0));

    fmt::print("...done\n");
}
