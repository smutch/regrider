#include <fmt/core.h>
#include <fmt/ostream.h>
#include <fstream>
#include <vector>
#include <string>

#include "../include/cxxopts.hpp"

static void read_gbptrees(const std::string fname_in) {
    fmt::print("Reading gbpTrees grid file {}...\n", fname_in);
    std::ifstream ifs(fname_in, std::ios::binary | std::ios::in);

    std::vector<int> n_cell(3);
    ifs.read(reinterpret_cast<char*>(n_cell.data()), sizeof(int) * 3);
    fmt::print("n_cell = {}\n", fmt::join(n_cell, ","));

    std::vector<double> box_size(3);
    ifs.read(reinterpret_cast<char*>(box_size.data()), sizeof(double) * 3);
    fmt::print("box_size = {}\n", fmt::join(box_size, ","));

    int32_t n_grids;
    ifs.read(reinterpret_cast<char*>(&n_grids), sizeof(int) * 3);
    fmt::print("n_cell = {}\n", n_grids);

    int32_t ma_scheme;
    ifs.read(reinterpret_cast<char*>(&ma_scheme), sizeof(int) * 3);
    fmt::print("n_cell = {}\n", ma_scheme);

    ifs.close();
    fmt::print("...done\n");
}

static void write_gbptrees(const std::string fname_out) {
    
}


int main(int argc, char *argv[])
{
    cxxopts::Options options("regrider", "Downsample gbpTrees and VELOCIraptor trees using FFTW");

    options.add_options()
        ("d,dim", "new grid dimension", cxxopts::value<uint32_t>())
        ("g,gbptrees", "input gbpTrees grid file", cxxopts::value<std::string>())
        ("v,velociraptor", "input VELOCIraptor grid file", cxxopts::value<std::string>());

    auto vm = options.parse(argc, argv);

    if (vm.count("gbptrees") && vm.count("velociraptor")) {
        fmt::print(stderr, "Must specify either gbpTrees or VELOCIraptor file. Not both...\n");
        return 1;
    }

    if (vm.count("gbptrees")) {
        read_gbptrees(vm["gbptrees"].as<std::string>());
    }

    return 0;
}
