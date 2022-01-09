#include "puppeteer.hpp"

#include "file_io/vtk_solution_writer.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_3d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_initial_conditions.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_parameters.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_solution_monitor.hpp"
#include "governing_equations/initial_conditions.hpp"
#include "governing_equations/residual_calculator.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "time_stepping/time_integrator.hpp"
#include "time_stepping/time_integrator_options.hpp"
#include "utils/registry.hpp"


Puppeteer::Puppeteer(const std::vector<std::string>& cmd_line)
{
  ProgramOptionsParser parser;
  parser.parseInputOptions(cmd_line);

  Logger::get().initLogFile(Options::get().log_file());
  Logger::get().InfoMessage(parser.getConfigurationFileTemplate());

  //
  // construct the primary simulation components
  //
  geom_               = geomFactory();
  rhs_                = rhsFactory(*geom_);
  initial_conditions_ = initialConditionsFactory();
  time_integrator_    = timeIntegratorFactory(*rhs_);

  //
  // construct and attach observers
  //
  attachResidualObservers(*time_integrator_, *rhs_);
  attachSolutionObservers(*time_integrator_, *geom_, *rhs_);
  attachTimeObservers(*time_integrator_);
}


Puppeteer::~Puppeteer() = default;


void
Puppeteer::run()
{
  profile();
  auto run_timer = Logger::get().timer("Starting solver");

  time_integrator_->solve(*rhs_, *initial_conditions_);

  const std::string vtk_output = [&]() {
    if (Options::get().solution_output_file().empty()) {
      return Options::get().log_file();
    }
    return Options::get().solution_output_file();
  }();

  write_solution_to_vtk(time_integrator_->getCurrentSolutionState(), *geom_, vtk_output);

  Logger::get().InfoMessage(Profiler::get().finalize());
}


std::unique_ptr<SpatialDiscretization>
Puppeteer::geomFactory()
{
  const int nhalo = 2;  // should be dependent on equation and order / stencil size

  const int dim = Options::get().dimension();
  switch (dim) {
    case 2:
      return std::make_unique<Discretization2DCart>(
          Options::get().domain_resolution(),
          nhalo,
          Options::get().domain_x_begin(),
          Options::get().domain_x_end(),
          Options::get().domain_x_begin());
    case 3:
      return std::make_unique<Discretization3DCart>(
          Options::get().domain_resolution(),
          nhalo,
          Options::get().domain_x_begin(),
          Options::get().domain_x_end(),
          Options::get().domain_x_begin(),
          Options::get().domain_x_begin());
    default:
      Logger::get().FatalMessage("Only dimensions 2 and 3 supported.");
  }
  return nullptr;
}


std::unique_ptr<TimeIntegrableRHS>
Puppeteer::rhsFactory(SpatialDiscretization& geom)
{
  Logger::get().FatalAssert(
      Options::get().equation_type() == "cahn-hilliard", "Only equation_type=cahn-hilliard currently supported.");

  CahnHilliardParameters params;

  const int dim = Options::get().dimension();
  switch (dim) {
    case 2:
      return std::make_unique<CahnHilliard2DFD>(dynamic_cast<Discretization2DCart&>(geom), params);
      ;
    case 3:
      return std::make_unique<CahnHilliard3DFD>(dynamic_cast<Discretization3DCart&>(geom), params);
      ;
    default:
      Logger::get().FatalMessage("Only dimensions 2 and 3 supported.");
  }

  return nullptr;
}


std::unique_ptr<InitialConditions>
Puppeteer::initialConditionsFactory()
{
  Logger::get().FatalAssert(
      Options::get().equation_type() == "cahn-hilliard", "Only equation_type=cahn-hilliard currently supported.");

  CahnHilliardParameters params;
  return std::make_unique<CahnHilliardInitialConditions>(params);
}


std::unique_ptr<TimeIntegrator>
Puppeteer::timeIntegratorFactory(TimeIntegrableRHS& rhs)
{
  auto time_opts = getTimeOptionsUsingGlobalOptions();
  return FactoryRegistry<TimeIntegrator>::get().lookup(Options::get().time_integrator())(rhs, time_opts);
}


void
Puppeteer::attachTimeObservers(TimeIntegrator& time_integrator)
{
  time_integrator.registerObserver(
      Event::TimeStepComplete, [&]() { Logger::get().setTimeMonitor(time_integrator.getIterationStatus()); });

  time_integrator.registerObserver(Event::TimeStepComplete, [&]() { Logger::get().updateLog(); });
}

void
Puppeteer::attachResidualObservers(TimeIntegrator& time_integrator, TimeIntegrableRHS& rhs)
{
  time_integrator.registerObserver(Event::TimeStepComplete, [&]() {
    Logger::get().setResidualMonitor(residual_calculator(rhs, time_integrator.getCurrentResidual()));
  });
}

void
Puppeteer::attachSolutionObservers(TimeIntegrator& time_integrator, SpatialDiscretization& geom, TimeIntegrableRHS& rhs)
{
  const int dim = Options::get().dimension();
  switch (dim) {
    case 2:
      time_integrator.registerObserver(Event::TimeStepComplete, [&]() {
        Logger::get().setSolutionMonitor(cahn_hilliard_solution_monitor(
            rhs, dynamic_cast<Discretization2DCart&>(geom), time_integrator.getCurrentSolutionState()));
      });
      return;
    case 3:
      time_integrator.registerObserver(Event::TimeStepComplete, [&]() {
        Logger::get().setSolutionMonitor(cahn_hilliard_solution_monitor(
            rhs, dynamic_cast<Discretization3DCart&>(geom), time_integrator.getCurrentSolutionState()));
      });
      return;
    default:
      Logger::get().FatalMessage("Only dimensions 2 and 3 supported.");
  }
}
