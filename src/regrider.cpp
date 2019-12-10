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
#include <cstring>
#include <fftw3.h>


static void read_gbptrees(const std::string fname_in, const std::string grid_name)
{
    fmt::print("Reading gbpTrees grid file {}...\n", fname_in);
    std::ifstream ifs(fname_in, std::ios::binary | std::ios::in);

    std::vector<int> n_cell(3);
    ifs.read((char*)(n_cell.data()), sizeof(int) * 3);
    fmt::print("n_cell = {}\n", fmt::join(n_cell, ","));

    std::vector<double> box_size(3);
    ifs.read((char*)(box_size.data()), sizeof(double) * 3);
    fmt::print("box_size = {}\n", fmt::join(box_size, ","));

    int32_t n_grids;
    ifs.read((char*)(&n_grids), sizeof(int));
    fmt::print("n_grids = {}\n", n_grids);

    int32_t ma_scheme;
    ifs.read((char*)(&ma_scheme), sizeof(int));
    fmt::print("ma_scheme = {}\n", ma_scheme);

    auto n_logical = n_cell[0] * n_cell[1] * n_cell[2];
    auto n_padded = n_cell[0] * n_cell[1] * 2*(n_cell[2]/2+1);
    std::unique_ptr<float, void(*)(float*)> orig(fftwf_alloc_real(n_padded), [](float* orig){ fftwf_free(orig); });

    bool found = false;
    for(int32_t ii=0; ii<n_grids && !found; ++ii){

        std::string ident(32, '\0');
        ifs.read((char*)(ident.data()), sizeof(ident));
        ident.resize(strlen(ident.c_str()));
        
        if (grid_name == ident) {
            fmt::print("Reading grid {}...\n", ident);
            ifs.read((char*)orig.get(), sizeof(float)*n_logical);
            found = true;
        } else {
            fmt::print("Skipping grid {}...\n", ident);
            ifs.seekg(sizeof(float)*n_logical, std::ifstream::cur);
        }

    }

    ifs.close();

    if (!found) {
        fmt::print(stderr, "Failed to find grid named `{}' in file!\n", grid_name);
    }

    // DEBUG
    {
        std::vector subset(orig.get(), orig.get() + 10);
        fmt::print("First 10 elements = {}\n", fmt::join(subset, ","));
    }

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

    if (vm.count("gbptrees")) {
        read_gbptrees(vm["gbptrees"].as<std::string>(), vm["name"].as<std::string>());
    }

    return 0;
}
