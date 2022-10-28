#pragma once

#include <vector>
#include <string>


class Options
{
  // only allow ProgramOptionsParser to populate the Options members
  friend class ProgramOptionsParser;

 public:
  // access to the singleton object
  static Options& get()
  {
    static Options instance;
    return instance;
  }

  // provide read-only access to all program options,
  // accessible only through the singleton.

  // file options
  const std::string& config_file() const { return config_file_; }
  const std::string& solution_output_file() const { return solution_output_file_; }
  const std::string& initial_condition_file() const { return initial_condition_file_; }
  const std::string& log_file() const { return log_file_; }
  int                log_frequency() const { return log_frequency_; }
  int                save_every() const { return save_every_; }

  // governing equations options
  const std::string& equation_type() const { return equation_type_; }
  double ch_m() const { return ch_m_; }
  double ch_eps2() const { return ch_eps2_; }
  double ch_sigma() const { return ch_sigma_; }
  bool   use_fused_kernels() const { return use_fused_kernels_; }


  // spatial_discretization options;
  double domain_x_begin() const { return domain_x_begin_; }
  double domain_x_end() const { return domain_x_end_; }
  int    domain_resolution() const { return domain_resolution_; }
  int    dimension() const { return dimension_; }
  const std::string& bc_x() const { return bc_x_; }
  const std::string& bc_y() const { return bc_y_; }
  const std::string& bc_z() const { return bc_z_; }
  
  // time stepping options
  const std::string& time_integrator() const { return time_integrator_; }
  double             time_step_size() const { return time_step_size_; }
  double             initial_time() const { return initial_time_; }
  double             final_time() const { return final_time_; }
  bool               use_cfl_time_step() const { return use_cfl_time_step_; }
  double             cfl() const { return cfl_; }
  bool               use_adaptive_time_step() const { return use_adaptive_time_step_; }
  double             min_time_step_size() const { return min_time_step_size_; }
  double             max_time_step_size() const { return max_time_step_size_; }
  unsigned int       max_time_steps() const { return max_time_steps_; }
  double             time_rel_err_tol() const { return time_rel_err_tol_; }
  double             converged_rel_tol() const { return converged_rel_tol_; }
  double             converged_abs_tol() const { return converged_abs_tol_; }

 private:
  Options() {}
  ~Options() {}

  // file options
  std::string config_file_;
  std::string solution_output_file_;
  std::string initial_condition_file_;
  std::string log_file_;
  int         log_frequency_;
  int         save_every_;

  // governing_equations options
  std::string         equation_type_;
  double              ch_m_;
  double              ch_eps2_;
  double              ch_sigma_;
  bool                use_fused_kernels_;

  // spatial_discretization options;
  double domain_x_begin_;
  double domain_x_end_;
  int domain_resolution_;
  int dimension_;
  std::string bc_x_;
  std::string bc_y_;
  std::string bc_z_;

  // time stepping options
  std::string time_integrator_;
  double      time_step_size_;
  double      initial_time_;
  double      final_time_;
  bool        use_cfl_time_step_;
  double      cfl_;
  bool        use_adaptive_time_step_;
  double      min_time_step_size_;
  double      max_time_step_size_;
  unsigned    max_time_steps_;
  double      time_rel_err_tol_;
  double      converged_rel_tol_;
  double      converged_abs_tol_;
};
