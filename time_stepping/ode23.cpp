#include "ode23.hpp"

ODE23::ODE23(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts)
    : TimeIntegrator(rhs, opts),
      controller_(opts.time_rel_err_tol_)
{
  k1_rhs_         = rhs.createSolutionState();
  k2_rhs_         = rhs.createSolutionState();
  k3_rhs_         = rhs.createSolutionState();
  k4_rhs_         = rhs.createSolutionState();
  stage_solution_ = rhs.createSolutionState();
  error_estimate_ = 1.0;
  first_step_     = true;
}


void
ODE23::doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
{
  profile();

  // first-same-as-last (FSAL) property.
  if (first_step_) {
    rhs.evalRHSImpl(state, time, *k1_rhs_);
    first_step_ = false;
  }
  else {
    for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
      auto k1 = write_access(k1_rhs_->getVec(ivec));
      auto k4 = read_access(k4_rhs_->getVec(ivec));

      maDGForAll(i, 0, k1.size(), {  //
        k1[i] = k4[i];
      });
    }
  }


  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto stage_sol = write_access(stage_solution_->getVec(ivec));
    auto sol0      = read_access(state.getVec(ivec));
    auto k1        = read_access(k1_rhs_->getVec(ivec));

    const real_wp a21 = a21_;

    maDGForAll(i, 0, stage_sol.size(), {  //
      stage_sol[i] = sol0[i] + dt * a21 * k1[i];
    });
  }


  rhs.evalRHSImpl(*stage_solution_, time + c2_ * dt, *k2_rhs_);


  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto stage_sol = write_access(stage_solution_->getVec(ivec));
    auto sol0      = read_access(state.getVec(ivec));
    auto k1        = read_access(k1_rhs_->getVec(ivec));
    auto k2        = read_access(k2_rhs_->getVec(ivec));

    const real_wp a31 = a31_;
    const real_wp a32 = a32_;

    maDGForAll(i, 0, stage_sol.size(), {  //
      stage_sol[i] = sol0[i] + dt * (a31 * k1[i] + a32 * k2[i]);
    });
  }


  rhs.evalRHSImpl(*stage_solution_, time + c3_ * dt, *k3_rhs_);


  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto residual = write_access(dstate_dt.getVec(ivec));
    auto k1       = read_access(k1_rhs_->getVec(ivec));
    auto k2       = read_access(k2_rhs_->getVec(ivec));
    auto k3       = read_access(k3_rhs_->getVec(ivec));

    const real_wp a41 = a41_;
    const real_wp a42 = a42_;
    const real_wp a43 = a43_;

    maDGForAll(i, 0, residual.size(), {  //
      residual[i] = a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
    });
  }


  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto sol = read_write_access(state.getVec(ivec));
    auto res = read_access(dstate_dt.getVec(ivec));

    maDGForAll(i, 0, sol.size(), {  //
      sol[i] += dt * res[i];
    });
  }


  // stage for error computation
  rhs.evalRHSImpl(state, time + c4_ * dt, *k4_rhs_);

  error_estimate_ = 0;
  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto k1  = read_access(k1_rhs_->getVec(ivec));
    auto k2  = read_access(k2_rhs_->getVec(ivec));
    auto k3  = read_access(k3_rhs_->getVec(ivec));
    auto k4  = read_access(k4_rhs_->getVec(ivec));
    auto res = read_access(dstate_dt.getVec(ivec));

    const real_wp d1 = d1_;
    const real_wp d2 = d2_;
    const real_wp d3 = d3_;
    const real_wp d4 = d4_;

    ReduceMaxReal max_residual(error_estimate_);

    maDGForAll(i, 0, res.size(), {  //
      const real_wp local_error_residual = d1 * k1[i] + d2 * k2[i] + d3 * k3[i] + d4 * k4[i];

      const real_wp local_error_estimate = dt * abs(res[i] - local_error_residual);

      max_residual.max(abs(local_error_estimate));
    });

    error_estimate_ = std::max(error_estimate_, max_residual.get());
  }
}


double
ODE23::computeNextDT(double current_dt)
{
  return controller_.adaptStep(error_estimate_, current_dt);
}


double
ODE23::PIDController::adaptStep(const double error_estimate, const double current_step)
{
  double err_n = std::max(error_estimate, std::numeric_limits<double>::epsilon());
  if (err_history.size() != 3) {
    step_history.push_back(current_step);
    err_history.push_back(err_n);
    return current_step;
  }
  else {
    err_history.pop_front();  // only keep 3 elements in the history deque
    step_history.pop_front();
    step_history.push_back(current_step);
    err_history.push_back(err_n);

    double wn    = step_history[2] / step_history[1];
    double alpha = (kI + kP + kD * (2.0 * wn) / (1.0 + wn)) / errest_order_;
    double beta  = (kP + 2.0 * wn * kD) / errest_order_;
    double gamma = kD * (2.0 * wn * wn) / (1.0 + wn) / errest_order_;

    double step = safety_ * current_step * pow(errtol_ / err_history[2], alpha) * pow(err_history[1] / errtol_, beta) *
                  pow(errtol_ / err_history[0], gamma);

    // restricts timestep from increasing by more than double or decreasing
    // by more than 1/10th between any two given timesteps.
    step = std::max(std::min(step, current_step * 1.1), current_step * 0.01);
    return step;
  }
}


static auto ode23instance = FactoryRegistry<TimeIntegrator>::get().add(
    "ode23",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) { return std::make_unique<ODE23>(rhs, opts); });
