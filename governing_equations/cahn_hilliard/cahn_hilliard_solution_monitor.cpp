#include "cahn_hilliard_solution_monitor.hpp"

#include "cahn_hilliard_state.hpp"
#include "data_structures/exec_includes.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "typedefs.hpp"

#include "utils/asymptotic_limit.hpp"



CahnHilliardSolutionMonitor::CahnHilliardSolutionMonitor(
    const TimeIntegrableRHS&     rhs,
    const SpatialDiscretization& geom)
    : rhs_(rhs), geom_(geom), iters_for_curve_fit_(100), iter_start_curve_fit_(1000)
{
  try {
    dynamic_cast<const Discretization2DCart&>(geom);
    is3D_ = false;
  }
  catch (...) {
    is3D_ = true;
  }
}


double CahnHilliardSolutionMonitor::asymptoticCurveFit(const std::deque<double>& time, const std::deque<double>& gradc, int iter)
{
  if (time.size() != iters_for_curve_fit_ || iter < iter_start_curve_fit_){
    return 1;
  }

  return asymptoticLimit(time, gradc);
}


std::vector<std::pair<std::string, double>>
CahnHilliardSolutionMonitor::eval(const SolutionState& state, double time, int iter)
{
  std::vector<std::pair<std::string, double>> monitors;

  if (is3D_) {
    const auto& g3d = dynamic_cast<const Discretization3DCart&>(geom_);
    monitors        = cahn_hilliard_solution_monitor_3d(rhs_, g3d, state);
  }
  else {
    const auto& g2d = dynamic_cast<const Discretization2DCart&>(geom_);
    monitors        = cahn_hilliard_solution_monitor_2d(rhs_, g2d, state);
  }


  time_history_.push_back(time);
  gradc_history_.push_back(monitors[2].second);

  if (time_history_.size() > iters_for_curve_fit_){ time_history_.pop_front(); }
  if (gradc_history_.size() > iters_for_curve_fit_){ gradc_history_.pop_front(); }

  monitors.push_back({"gradc_lim", asymptoticCurveFit(time_history_, gradc_history_, iter)});

  return monitors;
}


std::vector<std::pair<std::string, double>>
CahnHilliardSolutionMonitor::cahn_hilliard_solution_monitor_2d(
    const TimeIntegrableRHS&    rhs,
    const Discretization2DCart& geom,
    const SolutionState&        state)
{
  ReduceMinReal c_min(std::numeric_limits<real_wp>::max());
  ReduceMaxReal c_max(std::numeric_limits<real_wp>::min());
  ReduceSumReal mag_gradc_integrated(0);

  auto f   = read_access(dynamic_cast<const CahnHilliardState&>(state).c());
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
CahnHilliardSolutionMonitor::cahn_hilliard_solution_monitor_3d(
    const TimeIntegrableRHS&    rhs,
    const Discretization3DCart& geom,
    const SolutionState&        state)
{
  ReduceMinReal c_min(std::numeric_limits<real_wp>::max());
  ReduceMaxReal c_max(std::numeric_limits<real_wp>::min());
  ReduceSumReal mag_gradc_integrated(0);

  auto f   = read_access(dynamic_cast<const CahnHilliardState3D&>(state).c());
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
