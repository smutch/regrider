/*
 * regrider: Downsample gbpTrees and VELOCIraptor grids using FFTW
 * Copyright © 2021 Simon Mutch
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

#include <H5Cpp.h>
#include <array>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <vector>

#include "utils.hpp"
#include "velociraptor.hpp"

enum property
{
  X_VELOCITY,
  Y_VELOCITY,
  Z_VELOCITY,
  DENSITY
};

void regrid_velociraptor(const std::string fname_in, const std::string fname_out, const int new_dim)
{
  fmt::print("Regridding VELOCIraptor file {}\n", fname_in);

  auto file_in = H5::H5File(fname_in, H5F_ACC_RDONLY);
  auto file_out = H5::H5File(fname_out, H5F_ACC_RDWR);

  int _dim = 0;
  std::array<int, 3> new_n_cell = { new_dim, new_dim, new_dim };
  {
    auto group_in = file_in.openGroup("/Parameters");
    auto attr = group_in.openAttribute("DensityGrids:grid_dim");
    std::string _data;
    attr.read(attr.getDataType(), _data);
    _dim = std::stoi(_data);
  }
  std::array<int, 3> n_cell = { _dim, _dim, _dim };
  std::array<double, 3> box_size = { 0, 0, 0 };
  {
    auto attr = file_in.openGroup("/Header").openAttribute("BoxSize");
    attr.read(attr.getDataType(), box_size.data());
  }

  fmt::print("n_cell = [{}] --> [{}]\n", fmt::join(n_cell, ", "), fmt::join(new_n_cell, ", "));
  fmt::print("box_size = {:.2f}\n", fmt::join(box_size, ", "));

  auto grid = Grid(n_cell, box_size);
  const double radius = (double)grid.box_size[0] / (double)new_dim * 0.5;

  auto group_out = file_out.createGroup("/PartType1").createGroup("/Grids");
  auto group_in = file_in.openGroup("/PartType1/Grids");

  for (int property = X_VELOCITY; property <= DENSITY; ++property) {
    std::string dset_name;
    switch (property) {
      case X_VELOCITY:
        dset_name = "Vx";
        break;
      case Y_VELOCITY:
        dset_name = "Vy";
        break;
      case Z_VELOCITY:
        dset_name = "Vz";
        break;
      case DENSITY:
        dset_name = "Density";
        break;
      default:
        fmt::print(stderr, "Unrecognised grid property!");
        break;
    }

    // We do this here as the Grid may have already been subsampled in a
    // previous iteration.
    grid.update_properties(n_cell);

    {
      fmt::print("Reading grid {}... ", dset_name);
      auto dset = group_in.openDataSet(dset_name);
      dset.read(grid.get(), dset.getDataType());
      print_done();
    }

    grid.filter(Grid::filter_type::real_top_hat, radius);
    grid.sample(new_n_cell);

    fmt::print("Writing subsampled grid {}... ", dset_name);
    std::array<hsize_t, 3> dims = { static_cast<unsigned long long>(n_cell[0]),
                                    static_cast<unsigned long long>(n_cell[1]),
                                    static_cast<unsigned long long>(n_cell[2]) };
    group_out.createDataSet(dset_name, H5::PredType::NATIVE_FLOAT, H5::DataSpace(3, dims.data()));

    print_done();
  }

  // Remember to update the grid dimensions
  group_out = file_out.openGroup("/Parameters");
  group_out.openAttribute("DensityGrids:grid_dim").write(H5::PredType::NATIVE_INT, new_n_cell.data());
  group_out.openAttribute("Snapshots:grid_dim").write(H5::PredType::NATIVE_INT, new_n_cell.data());
}
