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

#ifndef VELOCIRAPTOR_H
#define VELOCIRAPTOR_H

#include "grid.hpp"
#include <string>

/** Regrid a VELOCIraptor file.
 *
 * @param fname_in The path to the input file to be regridded
 * @param fname_out The path to the new output file to be created
 * @param new_dim The new size of the grid (assuming cubic dimensions)
 */
void regrid_velociraptor(const std::string fname_in, const std::string fname_out, const int new_dim);

#endif
