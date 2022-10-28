#include "program_options_parser.hpp"

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <vector>
#include <iomanip>

#include "logger/logger.hpp"
#include "program_options.hpp"
#include "time_stepping/time_integrator.hpp"
#include "utils/registry.hpp"

#include <boost/algorithm/string/predicate.hpp>


namespace po = boost::program_options;

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
  for (size_t i = 0; i < vec.size(); ++i) {
    os << vec[i] << " ";
  }
  return os;
}

namespace std
{
std::ostream& operator<<(std::ostream& os, const std::vector<double>& vec)
{
  for (auto i : vec) {
    os << i << " ";
  }
  return os;
}
}


namespace std
{
std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vec)
{
  for (auto i : vec) {
    os << i << " ";
  }
  return os;
}
}


std::istream& operator>>(std::istream& in, LogLevel& level)
{
  std::string token;
  in >> token;
  if (boost::iequals(token, "trace")) {
    level = LogLevel::trace;
  } else if (boost::iequals(token, "debug")) {
    level = LogLevel::debug;
  } else if (boost::iequals(token, "info")) {
    level = LogLevel::info;
  } else if (boost::iequals(token, "warning")) {
    level = LogLevel::warning;
  } else if (boost::iequals(token, "error")) {
    level = LogLevel::error;
  } else if (boost::iequals(token, "fatal")) {
    level = LogLevel::fatal;
  } else {
    in.setstate(std::ios_base::failbit);
  }
  return in;
}


ProgramOptionsParser::ProgramOptionsParser()
{
  base.add_options()("help", "Produces this help message.");

  file.add_options()(
      "config_file",
      po::value<std::string>(&Options::get().config_file_),
      "Input configuration options file. Any command line options will override config file "
      "settings.")(
      "solution_output_file",
      po::value<std::string>(&Options::get().solution_output_file_),
      "Solution output file name, without file extension.")(
      "initial_condition_file",
      po::value<std::string>(&Options::get().initial_condition_file_)->default_value("none"),
      "Filename containing and initial solution. Overrides input values if not 'none'.")(
      "log_file",
      po::value<std::string>(&Options::get().log_file_)->default_value("maDGiCart"),
      "Solution output file name, without file extension.")(
      "log_frequency",
      po::value<int>(&Options::get().log_frequency_)->default_value(20),
      "Logger update frequency.")(
      "save_every",
      po::value<int>(&Options::get().save_every_)->default_value(0),
      "If not zero, save a solution every N time step.");

  physics.add_options()(
      "equation_type",
      po::value<std::string>(&Options::get().equation_type_)->default_value("cahn-hilliard"),
      "Governing equations to solve.")(
      "m",
      po::value<double>(&Options::get().ch_m_)->default_value(0.0),
      "Macroscopic average volume fraction of concentration.")(
      "eps2",
      po::value<double>(&Options::get().ch_eps2_)->default_value(2.419753086419753e-05),
      "epsilon^2, the biharmonic coefficient.")(
      "sigma",
      po::value<double>(&Options::get().ch_sigma_)->default_value(291.4488109293305),
      "sigma, the coefficient of the linearly stabilizing term about m.")(
      "use_fused_kernels",
      po::value<bool>(&Options::get().use_fused_kernels_)->default_value(false),
      "Use fused kernels for 3D CH");

  discretization.add_options()(
      "domain_x_begin",
      po::value<double>(&Options::get().domain_x_begin_)->default_value(0),
      "Minimum x value for Cartesian domain")(
      "domain_x_end",
      po::value<double>(&Options::get().domain_x_end_)->default_value(1),
      "Maximum x value for Cartesian domain")(
      "domain_resolution",
      po::value<int>(&Options::get().domain_resolution_)->default_value(128),
      "Cartesian mesh with [domain_resolution]^N nodes")(
      "dimension",
      po::value<int>(&Options::get().dimension_)->default_value(2),
      "Spatial dimension, 2 or 3.")(
      "bc_x",
      po::value<std::string>(&Options::get().bc_x_)->default_value("periodic"),
      "Boundary condition in x-direction, periodic or neumann. Only for 3D.")(
      "bc_y",
      po::value<std::string>(&Options::get().bc_y_)->default_value("periodic"),
      "Boundary condition in y-direction, periodic or neumann. Only for 3D.")(
      "bc_z",
      po::value<std::string>(&Options::get().bc_z_)->default_value("periodic"),
      "Boundary condition in z-direction, periodic or neumann. Only for 3D.");

  time_stepping.add_options()(
      "time_integrator",
      po::value<std::string>(&Options::get().time_integrator_)->default_value("rk3ssp"),
      "Time integration method.")(
      "time_step_size",
      po::value<double>(&Options::get().time_step_size_)->default_value(1.e-6),
      "Unsteady time step size. For adaptive methods this is the initial step size.")(
      "initial_time",
      po::value<double>(&Options::get().initial_time_)->default_value(0),
      "Initial time for unsteady.")(
      "final_time",
      po::value<double>(&Options::get().final_time_)->default_value(1),
      "Stopping time for unsteady.")(
      "max_time_steps",
      po::value<unsigned>(&Options::get().max_time_steps_)->default_value(10000),
      "Maximum number of unsteady time steps.")(
      "use_cfl_time_step",
      po::value<bool>(&Options::get().use_cfl_time_step_)->default_value(false),
      "Use CFL-based time stepping instead of dimensional time step.")(
      "cfl",
      po::value<double>(&Options::get().cfl_)->default_value(1.0),
      "CFL number for CFL-based time stepping.")(
      "use_adaptive_time_step",
      po::value<bool>(&Options::get().use_adaptive_time_step_)->default_value(false),
      "Use error-adaptive time stepping.")(
      "min_time_step_size",
      po::value<double>(&Options::get().min_time_step_size_)
          ->default_value(std::numeric_limits<double>::epsilon()*100),
      "For adaptive time stepping, the minimum time step size.")(
      "max_time_step_size",
      po::value<double>(&Options::get().max_time_step_size_)
          ->default_value(std::numeric_limits<double>::max()/100),
      "For adaptive time stepping, the maximum time step size.")(
      "time_rel_err_tol",
      po::value<double>(&Options::get().time_rel_err_tol_)->default_value(1e-2),
      "Relative error tolerance for adaptive time stepping.")(
      "converged_rel_tol",
      po::value<double>(&Options::get().converged_rel_tol_)->default_value(0),
      "Relative residual to consider solution converged. 0 indicates use max time or max timesteps.")(
      "converged_abs_tol",
      po::value<double>(&Options::get().converged_abs_tol_)->default_value(0),
      "Absolute residual to consider solution converged. 0 indicates use max time or max timesteps.");

}

void ProgramOptionsParser::parseInputOptions(const std::vector<std::string>& cmd_line)
{
  vm.clear();

  base.add(file);
  base.add(discretization);
  base.add(physics);
  base.add(time_stepping);
  visible.add(base);


  po::positional_options_description p;
  p.add("config_file", -1);
  po::store(
      po::command_line_parser(cmd_line)
          .options(base)
          .style(po::command_line_style::unix_style ^ po::command_line_style::allow_short)
          .positional(p)
          .run(),
      vm);

  po::notify(vm);
  std::ifstream ifs(Options::get().config_file_.c_str());
  po::store(po::parse_config_file(ifs, base), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << visible << "\n";
    printRegistryContents(std::cout);
    std::stringstream ss;
    std::cout << getConfigurationFileTemplate();
    exit(0);
    return;
  }
}


std::ostream& operator<<(std::ostream& os, const boost::any& a)
{
  // must determine the type of the boost::any value stored by program options
  // in order to use the << operator.
  using boost::any_cast;
  try {
    std::string str = any_cast<std::string>(a);
    os << str;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    double dbl = any_cast<double>(a);
    os << dbl;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    int in = any_cast<int>(a);
    os << in;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    bool b = any_cast<bool>(a);
    os << b;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    unsigned uin = any_cast<unsigned>(a);
    os << uin;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    std::vector<double> vec = any_cast<std::vector<double> >(a);
    os << vec;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  try {
    std::vector<std::string> vec = any_cast<std::vector<std::string> >(a);
    os << vec;
    return os;
  } catch (boost::bad_any_cast&) {
  }
  return os;
}

void ProgramOptionsParser::printRegistryContents(std::ostream& out) const
{
  out << FactoryRegistry<TimeIntegrator>::get() << std::endl;
}

std::string ProgramOptionsParser::getConfigurationFileTemplate() const
{
  std::stringstream out;

  out.setf(std::ios_base::boolalpha);

  out << "\n# ****** File Options ****** #\n\n";
  for (const auto& i : file.options()) {
    std::string varname = i->long_name();
    out << std::setw(30) << varname << " = " << vm[varname].value() << std::endl;
  }

  out << "\n# ****** Physics Options ****** #\n\n";
  for (const auto& i : physics.options()) {
    std::string varname = i->long_name();
    out << std::setw(30) << varname << " = " << vm[varname].value() << std::endl;
  }

  out << "\n# ****** Discretization Options ****** #\n\n";
  for (const auto& i : discretization.options()) {
    std::string varname = i->long_name();
    out << std::setw(30) << varname << " = " << vm[varname].value() << std::endl;
  }

  out << "\n# ****** Time Stepping Options ****** #\n\n";
  for (const auto& i : time_stepping.options()) {
    std::string varname = i->long_name();
    out << std::setw(30) << varname << " = " << vm[varname].value() << std::endl;
  }

  out << "\n";
  return out.str();
}
