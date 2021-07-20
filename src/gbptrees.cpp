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
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fstream>
#include <vector>

#include "gbptrees.hpp"
#include "utils.hpp"

void regrid_gbptrees(const std::string fname_in, const std::string fname_out, const int new_dim)
{
  fmt::print("Regridding gbpTrees file {}\n", fname_in);
  std::ifstream ifs(fname_in, std::ios::binary | std::ios::in);
  std::ofstream ofs(fname_out, std::ios::binary | std::ios::out);

  std::array<int, 3> n_cell;
  std::array<int, 3> new_n_cell = { new_dim, new_dim, new_dim };
  ifs.read((char*)(n_cell.data()), sizeof(int) * 3);
  fmt::print("n_cell = [{}] --> [{}]\n", fmt::join(n_cell, ", "), fmt::join(new_n_cell, ", "));
  ofs.write((char*)(new_n_cell.data()), sizeof(int) * 3);

  std::array<double, 3> box_size;
  ifs.read((char*)(box_size.data()), sizeof(double) * 3);
  fmt::print("box_size = {}\n", fmt::join(box_size, ","));
  ofs.write((char*)(box_size.data()), sizeof(double) * 3);

  int32_t n_grids;
  ifs.read((char*)(&n_grids), sizeof(int));
  fmt::print("n_grids = {}\n", n_grids);
  ofs.write((char*)(&n_grids), sizeof(int));

  int32_t ma_scheme;
  ifs.read((char*)(&ma_scheme), sizeof(int));
  fmt::print("ma_scheme = {}\n", ma_scheme);
  ofs.write((char*)(&ma_scheme), sizeof(int));

  auto grid = Grid(n_cell, box_size);
  const double radius = (double)grid.box_size[0] / (double)new_dim * 0.5;

  for (int ii = 0; ii < n_grids; ++ii) {

    std::string ident(32, '\0');
    ifs.read((char*)(ident.data()), sizeof(ident));
    ofs.write((char*)(ident.data()), sizeof(ident));

    ident.resize(strlen(ident.c_str()));
    fmt::print("\nGrid {}\n=================\n", ident);

    // We do this here as the Grid may have already been subsampled in a
    // previous iteration.
    grid.update_properties(n_cell);

    fmt::print("Reading grid... ");
    ifs.read((char*)grid.get(), sizeof(float) * grid.n_logical);
    print_done();

#ifdef DEBUG
    {
      std::vector<float> subset(grid.get(), grid.get() + 10);
      fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }
#endif

    grid.filter(Grid::filter_type::real_top_hat, radius);

#ifdef DEBUG
    {
      std::vector<float> subset(grid.get(), grid.get() + 10);
      fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }
#endif

    grid.sample(new_n_cell);

#ifdef DEBUG
    {
      std::vector<float> subset(grid.get(), grid.get() + 10);
      fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }
#endif

    fmt::print("Writing subsampled grid... ");
    ofs.write((char*)grid.get(), sizeof(float) * grid.n_logical);
    print_done();
  }

  ofs.close();
  ifs.close();

  print_done();
}
