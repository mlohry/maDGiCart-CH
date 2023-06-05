#pragma once

#include "governing_equations/time_integrable_rhs.hpp"

#include <deque>

class ODE23Smoother {
 public:
  ODE23Smoother(TimeIntegrableRHS& rhs);

  void doSmooth(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, real_wp time);

  real_wp stepSize() const { return dt_; }

  const SolutionState& getLastResidual() const { return *k4_rhs_; }

 private:
  real_wp dt_ = 1e-6;

  std::unique_ptr<SolutionState> k1_rhs_;
  std::unique_ptr<SolutionState> k2_rhs_;
  std::unique_ptr<SolutionState> k3_rhs_;
  std::unique_ptr<SolutionState> k4_rhs_;
  std::unique_ptr<SolutionState> stage_solution_;

  real_wp error_estimate_;

  // Butcher tableau coefficients
  const real_wp c2_ = 1. / 2.;
  const real_wp c3_ = 3. / 4.;
  const real_wp c4_ = 1.;

  const real_wp a21_ = 1. / 2.;

  const real_wp a31_ = 0.;
  const real_wp a32_ = 3. / 4.;

  const real_wp a41_ = 2. / 9.;
  const real_wp a42_ = 1. / 3.;
  const real_wp a43_ = 4. / 9.;

  //  const real_wp b1_ = 2. / 9.;
  //  const real_wp b2_ = 1. / 3.;
  //  const real_wp b3_ = 4. / 9.;
  //  const real_wp b4_ = 0;

  const real_wp d1_ = 7. / 24.;
  const real_wp d2_ = 1. / 4.;
  const real_wp d3_ = 1. / 3.;
  const real_wp d4_ = 1. / 8.;


  double cfl_scale_;

  real_wp computeNextDT(real_wp current_dt);


  class PIDController {
   public:
    PIDController()
    {
    }

    real_wp adaptStep(const real_wp error_estimate, const real_wp current_step);

   private:
    const real_wp errtol_ = 1.e-6;
    const real_wp errest_order_ = 2;
    const real_wp safety_       = 1.0;
    /*
     * Controller gains for PID timestep controller taken from
     * "Additive Runge-Kutta Schemes for Convection-Diffusion-Reaction Equations"
     * Kennedy & Carpenter, NASA/TM-2001-211038
     */
    const real_wp kI = 0.25, kP = 0.14, kD = 0.10;

    const double minstep_ = 1.e-8;
    const double maxstep_ = 100;


    std::deque<real_wp> step_history, err_history;
  } controller_;
};
