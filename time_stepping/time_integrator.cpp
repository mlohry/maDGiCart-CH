#include "time_integrator.hpp"

#include "governing_equations/residual_calculator.hpp"

TimeIntegrator::TimeIntegrator(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts)
    : rhs_(rhs),
      time_options_(opts),
      solution_state_(rhs.createSolutionState()),
      residual_state_(rhs.createSolutionState())
{
}


void TimeIntegrator::solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions)
{
  profile();
  auto& state = *solution_state_;

  initial_conditions.set(rhs, state);

  Logger::get().FatalAssert(state.nvecs() == rhs.nEquations(), "state.nvecs() == device_rhs.nEquations()");

  std::tie(time_, dt_, cfl_, iter_) = initialTimestep(rhs, time_options_, state);

  bool last_step        = false;
  bool max_time_reached = false;

  while (!max_time_reached) {

    doTimeStep(rhs, state, *residual_state_, time_, dt_);

    time_ += dt_;
    iter_++;

    notifyObservers(Event::TimeStepComplete);
    notifyObservers(Event::SolutionUpdate);

    if (last_step) {
      break;
    }

    std::tie(dt_, last_step, max_time_reached) =
        computeNextTimestep(rhs, time_options_, state, cfl_, time_, dt_, iter_);
  }
}


std::tuple<real_t, real_t, real_t, int_t>
TimeIntegrator::initialTimestep(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& options, const SolutionState& solution)
{
  real_t time      = options.t0_;
  real_t cfl       = options.cfl_;
  int_t  iteration = options.initial_step_;
  real_t dt = -1;

  if (options.use_cfl_time_step_) {
    Logger::get().FatalMessage("use_cfl_time_step CFL time stepping not implemented.");
  }
  else {
    dt = options.dt_initial_;
  }

  return std::make_tuple(time, dt, cfl, iteration);
}


double TimeIntegrator::computeNextDT(double current_dt)
{
  return current_dt;
}


std::tuple<real_t, bool, bool>
TimeIntegrator::computeNextTimestep(
    TimeIntegrableRHS&           rhs,
    const TimeIntegratorOptions& options,
    const SolutionState&         solution,
    real_t                       cfl,
    real_t                       current_time,
    real_t                       current_dt,
    int_t                        current_iter)
{
  bool   last_iteration   = false;
  bool   max_time_reached = false;
  real_t dt               = computeNextDT(current_dt);

  if (options.use_cfl_time_step_) {
    Logger::get().FatalMessage("use_cfl_time_step CFL time stepping not implemented.");
  }

  if (current_time + dt >= options.tfinal_) {
    dt             = options.tfinal_ - current_time;
    last_iteration = true;
  }

  // max_time_steps_=0 means no limit on iteration count
  if (options.max_time_steps_ && current_iter >= options.max_time_steps_) {
    max_time_reached = true;
  }

  if (options.converged_rel_tol_){
    auto resmon = Logger::get().getResidualMonitor();
    for (const auto& kv : resmon){
      if (kv.first == "rel_res"){
        if (kv.second < options.converged_rel_tol_)
        {
          Logger::get().InfoMessage("Exiting due to converged_rel_tol criterion.");
          max_time_reached = true;
        }
        break;
      }
    }
  }

  if (options.converged_abs_tol_){
    auto resmon = Logger::get().getResidualMonitor();
    if (resmon.front().second < options.converged_abs_tol_){
      Logger::get().InfoMessage("Exiting due to converged_abs_tol criterion.");
      max_time_reached = true;
    }
  }

  return std::make_tuple(dt, last_iteration, max_time_reached);
}