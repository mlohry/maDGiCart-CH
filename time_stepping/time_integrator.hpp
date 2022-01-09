#pragma once

#include "governing_equations/initial_conditions.hpp"
#include "governing_equations/time_integrable_rhs.hpp"
#include "time_integrator_options.hpp"
#include "typedefs.hpp"
#include "utils/observer.hpp"
#include "utils/registry.hpp"


class TimeIntegrator : public Observable {
 public:
  using CtorType = std::function<std::unique_ptr<TimeIntegrator>(TimeIntegrableRHS&, const TimeIntegratorOptions&)>;

  TimeIntegrator(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);

  virtual ~TimeIntegrator() {}

  void solve(TimeIntegrableRHS&, InitialConditions&);

  const SolutionState& getCurrentResidual() const { return *residual_state_; }
  const SolutionState& getCurrentSolutionState() const { return *solution_state_; }


  virtual std::vector<std::pair<std::string, double>> getIterationStatus() const
  {
    std::vector<std::pair<std::string, double>> status;

    status.push_back({"iter", getCurrentStep()});
    status.push_back({"time", getCurrentTime()});
    status.push_back({"dt", getTimeStepSize()});

    return status;
  }

  virtual void doTimeStep(TimeIntegrableRHS&, SolutionState&, SolutionState&, double time, double dt) = 0;

 protected:

  virtual int_t  getCurrentStep() const { return iter_; }
  virtual real_t getCurrentTime() const { return time_; }
  virtual real_t getTimeStepSize() const { return dt_; }

 private:
  TimeIntegrableRHS& rhs_;

  const TimeIntegratorOptions time_options_;
  int_t                       iter_;
  real_t                      time_;
  real_t                      dt_;
  real_t                      cfl_;

  std::unique_ptr<SolutionState> solution_state_;
  std::unique_ptr<SolutionState> residual_state_;


  /// returns (time0, dt0, cfl, iteration)
  std::tuple<real_t, real_t, real_t, int_t>
  initialTimestep(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& options, const SolutionState& solution);

  /// return (dt, do one more step then quit, quit now)
  std::tuple<real_t, bool, bool> computeNextTimestep(
      TimeIntegrableRHS&           rhs,
      const TimeIntegratorOptions& options,
      const SolutionState&         solution,
      real_t                       cfl,
      real_t                       current_time,
      real_t                       current_dt,
      int_t                        current_iter);
};
