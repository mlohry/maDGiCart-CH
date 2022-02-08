#pragma once

#include "time_integrator.hpp"

#include <deque>

/**
 * Bogacki-Shampine explicit Runge-Kutta method, aka ode23 in Matlab.
 *
 * 3rd order with embedded 2nd order error estimate, 4-stages.
 * Butcher tableau has FSAL property so it uses 3 function evaluations per
 * time step.
 */
class ODE23 : public TimeIntegrator {
 public:
  ODE23(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);

  void doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
      override;


  std::vector<std::pair<std::string, double>> getIterationStatus() const override
  {
    std::vector<std::pair<std::string, double>> status;

    status.push_back({"iter", getCurrentStep()});
    status.push_back({"time", getCurrentTime()});
    status.push_back({"dt", getTimeStepSize()});
    status.push_back({"time_err", error_estimate_});

    return status;
  }

 private:
  std::unique_ptr<SolutionState> k1_rhs_;
  std::unique_ptr<SolutionState> k2_rhs_;
  std::unique_ptr<SolutionState> k3_rhs_;
  std::unique_ptr<SolutionState> k4_rhs_;
  std::unique_ptr<SolutionState> stage_solution_;

  bool   first_step_;
  double error_estimate_;

  // Butcher tableau coefficients
  const double c2_ = 1. / 2.;
  const double c3_ = 3. / 4.;
  const double c4_ = 1.;

  const double a21_ = 1. / 2.;

  const double a31_ = 0.;
  const double a32_ = 3. / 4.;

  const double a41_ = 2. / 9.;
  const double a42_ = 1. / 3.;
  const double a43_ = 4. / 9.;

  const double b1_ = 2. / 9.;
  const double b2_ = 1. / 3.;
  const double b3_ = 4. / 9.;
  const double b4_ = 0;

  const double d1_ = 7. / 24.;
  const double d2_ = 1. / 4.;
  const double d3_ = 1. / 3.;
  const double d4_ = 1. / 8.;


  double computeNextDT(double current_dt) override;


  class PIDController {
   public:
    PIDController(TimeIntegratorOptions opts)
        : errtol_(opts.time_rel_err_tol_), minstep_(opts.min_time_step_size_), maxstep_(opts.max_time_step_size_)
    {
    }

    double adaptStep(const double error_estimate, const double current_step);

   private:
    const double errtol_;
    const double minstep_;
    const double maxstep_;
    const double errest_order_ = 2;
    const double safety_       = 1.0;
    /*
     * Controller gains for PID timestep controller taken from
     * "Additive Runge-Kutta Schemes for Convection-Diffusion-Reaction Equations"
     * Kennedy & Carpenter, NASA/TM-2001-211038
     */
    const double kI = 0.25, kP = 0.14, kD = 0.10;


    std::deque<double> step_history, err_history;
  } controller_;
};
