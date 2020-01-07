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

#include <cxxopts.hpp>
#include <fftw3.h>
#include <fmt/core.h>
#include <fstream>

#include "gbptrees.hpp"

int main(int argc, char* argv[])
{
    cxxopts::Options options(
        "regrider", "Downsample gbpTrees and VELOCIraptor trees using FFTW");

    options.add_options() // clang-format off
        ("d,dim", "new grid dimension", cxxopts::value<int>())
        ("g,gbptrees", "input gbpTrees grid file", cxxopts::value<std::string>())
        ("v,velociraptor", "input VELOCIraptor grid file", cxxopts::value<std::string>())
        ("o,output", "output file name", cxxopts::value<std::string>())
        ("h,help", "show help", cxxopts::value<bool>());

    auto vm = options.parse(argc, argv);

    if (vm.count("help")) {
        fmt::print(options.help());
    }

    if (vm.count("gbptrees") && vm.count("velociraptor")) {
        fmt::print(stderr, "Must specify either gbpTrees or VELOCIraptor file. Not both...\n");
        return 1;
    }

    if (!vm.count("dim")) {
        fmt::print(stderr, "Must specify new grid dimension...\n");
        return 1;
    }

    fftwf_init_threads();

    if (vm.count("gbptrees")) {
        regrid_gbptrees(vm["gbptrees"].as<std::string>(), vm["output"].as<std::string>(), vm["dim"].as<int>());
    }

    fftwf_cleanup_threads();

    return 0;
}
