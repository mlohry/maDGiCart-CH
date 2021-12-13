#pragma once

#include "time_integrator.hpp"

/**
 * 3-stage 3rd order explicit Runge-Kutta.
 *
 * "Strong Stability-Preserving High-Order Time Discretization Methods"
 * Gottlieb, Shu, and Tadmor, SIAM Review 2001 Vol. 43,No. 1,pp. 89–112.
 */
class RK3SSP : public TimeIntegrator {
 public:
  RK3SSP(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);

  void doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
      override;


 private:
  std::unique_ptr<SolutionState> k1_rhs_;
  std::unique_ptr<SolutionState> k2_rhs_;
  std::unique_ptr<SolutionState> k3_rhs_;
  std::unique_ptr<SolutionState> stage_solution_;


  // Butcher tableau coefficients
  const double c2_ = 1.0;
  const double c3_ = 0.5;

  const double a21_ = 1.0;
  const double a31_ = 0.25;
  const double a32_ = 0.25;
  const double b1_  = 1.0 / 6.0;
  const double b2_  = 1.0 / 6.0;
  const double b3_  = 2.0 / 3.0;
};
