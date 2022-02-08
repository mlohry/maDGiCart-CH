#pragma once

#include <string>
#include "typedefs.hpp"
#include "program_options/program_options.hpp"

// todo remove this? bad idea to have multiply defined variables
class TimeIntegratorOptions
{
 public:
  TimeIntegratorOptions()
  {
  }

  real_t      t0_ ;
  real_t      tfinal_;
  real_t      dt_initial_;
  bool        use_cfl_time_step_;
  int_t       max_time_steps_;
  int_t       initial_step_ = 0;
  real_t      cfl_;
  real_t      time_rel_err_tol_;
  bool        use_adaptive_time_step_;
  real_t      min_time_step_size_;
  real_t      max_time_step_size_;
  real_t      converged_rel_tol_;
  real_t      converged_abs_tol_;

  /**
   * Named parameter idiom. Construct TimeIntegratorOptions like
   *
   * TimeIntegratorOptions.t0(0).tfinal(1).dt_initial(0.1).filter_after_time_step(true)....
   */
  TimeIntegratorOptions& t0(real_t t)
  {
    t0_ = t;
    return *this;
  }
  TimeIntegratorOptions& tfinal(real_t t)
  {
    tfinal_ = t;
    return *this;
  }
  TimeIntegratorOptions& dt_initial(real_t dt)
  {
    dt_initial_ = dt;
    return *this;
  }
  TimeIntegratorOptions& use_cfl_time_step(bool b)
  {
    use_cfl_time_step_ = b;
    return *this;
  }
  TimeIntegratorOptions& max_time_steps(int_t maxts)
  {
    max_time_steps_ = maxts;
    return *this;
  }
  TimeIntegratorOptions& cfl(real_t c)
  {
    cfl_ = c;
    return *this;
  }
  TimeIntegratorOptions& rel_err_tol(real_t r)
  {
    time_rel_err_tol_ = r;
    return *this;
  }
  TimeIntegratorOptions& use_adaptive_time_step(bool b)
  {
    use_adaptive_time_step_ = b;
    return *this;
  }
  TimeIntegratorOptions& min_time_step_size(real_t dt)
  {
    min_time_step_size_ = dt;
    return *this;
  }
  TimeIntegratorOptions& max_time_step_size(real_t dt)
  {
    max_time_step_size_ = dt;
    return *this;
  }
  TimeIntegratorOptions& converged_rel_tol(real_t r)
  {
    converged_rel_tol_ = r;
    return *this;
  }
  TimeIntegratorOptions& converged_abs_tol(real_t a)
  {
    converged_abs_tol_ = a;
    return *this;
  }
};


inline TimeIntegratorOptions getTimeOptionsUsingGlobalOptions()
{
  return TimeIntegratorOptions()
      .t0(Options::get().initial_time())
      .tfinal(Options::get().final_time())
      .dt_initial(Options::get().time_step_size())
      .use_cfl_time_step(Options::get().use_cfl_time_step())
      .max_time_steps(Options::get().max_time_steps())
      .cfl(Options::get().cfl())
      .rel_err_tol(Options::get().time_rel_err_tol())
      .use_adaptive_time_step(Options::get().use_adaptive_time_step())
      .min_time_step_size(Options::get().min_time_step_size())
      .max_time_step_size(Options::get().max_time_step_size())
      .converged_rel_tol(Options::get().converged_rel_tol())
      .converged_abs_tol(Options::get().converged_abs_tol());
}
