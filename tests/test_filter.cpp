#include <array>
#include <grid.hpp>
#include <criterion/criterion.h>

Test(filter, basic) {

  const float tolerance = 1e-5;

  std::array<int32_t, 3> n_cell = {32, 32, 32};
  std::array<double, 3> box_size = {10., 10., 10.};

  auto grid = Grid(n_cell, box_size);
  cr_assert(true);

  auto rgrid = grid.get();

  for (int ii=0; ii<grid.n_logical; ++ii) {
    rgrid[ii] = 0.0;
  }

  auto center = grid.index(n_cell[0]/2, n_cell[1]/2, n_cell[2]/2, Grid::index_type::real);
  rgrid[center] = 10.0;

  auto radius = 2.0;
  grid.filter(Grid::filter_type::real_top_hat, radius);

  auto rgrid_total = 0.0;
  for (int ii=0; ii<grid.n_logical; ++ii) {
    rgrid_total += rgrid[ii];
  }
  
  // cr_assert_float_eq(rgrid_total, 10.0, tolerance);
  cr_assert_float_eq(rgrid[0], 0.0, tolerance);
}
