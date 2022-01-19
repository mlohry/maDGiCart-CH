#pragma once

#include "governing_equations/solution_monitor.hpp"

#include <deque>
#include <string>
#include <utility>
#include <vector>

class TimeIntegrableRHS;
class SolutionState;
class Discretization2DCart;
class Discretization3DCart;
class SpatialDiscretization;


class CahnHilliardSolutionMonitor : public SolutionMonitor {
 public:
  CahnHilliardSolutionMonitor(const TimeIntegrableRHS& rhs, const SpatialDiscretization& geom);

  std::vector<std::pair<std::string, double>> eval(const SolutionState& state, double time, int iter) override;

  std::vector<std::pair<std::string, double>> cahn_hilliard_solution_monitor_2d(
      const TimeIntegrableRHS&    rhs,
      const Discretization2DCart& geom,
      const SolutionState&        state);

  std::vector<std::pair<std::string, double>> cahn_hilliard_solution_monitor_3d(
      const TimeIntegrableRHS&    rhs,
      const Discretization3DCart& geom,
      const SolutionState&        state);

 private:
  const TimeIntegrableRHS&     rhs_;
  const SpatialDiscretization& geom_;
  const size_t                 iters_for_curve_fit_;
  const int                    iter_start_curve_fit_;
  bool                         is3D_;

  std::deque<double> time_history_;
  std::deque<double> gradc_history_;

  double asymptoticCurveFit(const std::deque<double>& time, const std::deque<double>& gradc, int iter);
};
