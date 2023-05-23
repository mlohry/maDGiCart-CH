#include "cahn_hilliard_solution_monitor.hpp"

#include "data_structures/exec_includes.hpp"
#include "data_structures/scalar_solution_state.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "typedefs.hpp"


std::vector<std::pair<std::string, double>>
cahn_hilliard_solution_monitor(
    const TimeIntegrableRHS&    rhs,
    const Discretization2DCart& geom,
    const SolutionState&        state)
{
  ReduceMinReal c_min(std::numeric_limits<real_wp>::max());
  ReduceMaxReal c_max(std::numeric_limits<real_wp>::min());
  ReduceSumReal mag_gradc_integrated(0);

  auto f   = read_access(dynamic_cast<const ScalarSolutionState2D&>(state).c());
  auto idx = read_access(geom.interiorIndices());

  const real_wp dx = geom.dx();
  const int     ni = geom.ni();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    f.getIJ(idx[ii], i, j);
    const real_wp cell_value = f(i, j);
    c_min.min(cell_value);
    c_max.max(cell_value);

    real_wp dfdx = 0;
    real_wp dfdy = 0;
    if (i == 0) {
      dfdx = (f(i + 1, j) - cell_value) / dx;
    }
    else if (i == ni - 1) {
      dfdx = (cell_value - f(i - 1, j)) / dx;
    }
    else {
      dfdx = (f(i + 1, j) - f(i - 1, j)) / (2.0 * dx);
    }

    if (j == 0) {
      dfdy = (f(i, j + 1) - cell_value) / dx;
    }
    else if (j == ni - 1) {
      dfdy = (cell_value - f(i, j - 1)) / dx;
    }
    else {
      dfdy = (f(i, j + 1) - f(i, j - 1)) / (2.0 * dx);
    }

    //  integral of || grad(c) || on the cell
    mag_gradc_integrated += sqrt(pow(dfdx, 2.0) + pow(dfdy, 2.0)) * dx * dx;
  });

  return std::vector<std::pair<std::string, double>>{
      {"min_c", c_min.get()}, {"max_c", c_max.get()}, {"grad_c_mag", mag_gradc_integrated.get()}};
}


std::vector<std::pair<std::string, double>>
cahn_hilliard_solution_monitor(
    const TimeIntegrableRHS&    rhs,
    const Discretization3DCart& geom,
    const SolutionState&        state)
{
  ReduceMinReal c_min(std::numeric_limits<real_wp>::max());
  ReduceMaxReal c_max(std::numeric_limits<real_wp>::min());
  ReduceSumReal mag_gradc_integrated(0);

  auto f   = read_access(dynamic_cast<const ScalarSolutionState3D&>(state).c());
  auto idx = read_access(geom.interiorIndices());

  const real_wp dx = geom.dx();
  const int     ni = geom.ni();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    int k;
    f.getIJK(idx[ii], i, j, k);
    const real_wp cell_value = f(i, j, k);
    c_min.min(cell_value);
    c_max.max(cell_value);

    real_wp dfdx = 0;
    real_wp dfdy = 0;
    real_wp dfdz = 0;

    if (i == 0) {
      dfdx = (f(i + 1, j, k) - cell_value) / dx;
    }
    else if (i == ni - 1) {
      dfdx = (cell_value - f(i - 1, j, k)) / dx;
    }
    else {
      dfdx = (f(i + 1, j, k) - f(i - 1, j, k)) / (2.0 * dx);
    }

    if (j == 0) {
      dfdy = (f(i, j + 1, k) - cell_value) / dx;
    }
    else if (j == ni - 1) {
      dfdy = (cell_value - f(i, j - 1, k)) / dx;
    }
    else {
      dfdy = (f(i, j + 1, k) - f(i, j - 1, k)) / (2.0 * dx);
    }

    if (k == 0) {
      dfdz = (f(i, j, k + 1) - cell_value) / dx;
    }
    else if (k == ni - 1) {
      dfdz = (cell_value - f(i, j, k - 1)) / dx;
    }
    else {
      dfdz = (f(i, j, k + 1) - f(i, j, k - 1)) / (2.0 * dx);
    }


    //  integral of || grad(c) || on the cell
    mag_gradc_integrated += sqrt(pow(dfdx, 2.0) + pow(dfdy, 2.0) + pow(dfdz, 2.0)) * dx * dx;
  });

  return std::vector<std::pair<std::string, double>>{
      {"min_c", c_min.get()}, {"max_c", c_max.get()}, {"grad_c_mag", mag_gradc_integrated.get()}};
}
