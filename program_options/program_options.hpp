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
  const std::string& log_file() const { return log_file_; }
  const std::string& solution_output_file() const { return solution_output_file_; }
  bool interactive_plots() const { return interactive_plots_; }

  // governing equations options
  const std::string& equation_type() const { return equation_type_; }


  // spatial_discretization options;
  double domain_x_begin() const { return domain_x_begin_; }
  double domain_x_end() const { return domain_x_end_; }
  int    domain_resolution() const { return domain_resolution_; }
  
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
  double             time_abs_err_tol() const { return time_abs_err_tol_; }
  double             time_rel_err_tol() const { return time_rel_err_tol_; }

 private:
  Options() {}
  ~Options() {}

  // file options
  std::string config_file_;
  std::string log_file_;
  std::string solution_output_file_;
  bool interactive_plots_;

  // governing_equations options
  std::string         equation_type_;

  // spatial_discretization options;
  double domain_x_begin_;
  double domain_x_end_;
  int domain_resolution_;

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
  double      time_abs_err_tol_;
  double      time_rel_err_tol_;
};
