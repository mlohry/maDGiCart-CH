#include "puppeteer.hpp"

#include <matplot/matplot.h>
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_initial_conditions.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_parameters.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_solution_monitor.hpp"
#include "governing_equations/initial_conditions.hpp"
#include "governing_equations/residual_calculator.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
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

  if (Options::get().interactive_plots()) {
    plotSolution(time_integrator_->getCurrentSolutionState());
  }
}

void
Puppeteer::plotSolution(const SolutionState& state)
{
  using vec2d = std::vector<std::vector<double>>;
  const int N = geom_->ni();

  vec2d X, Y, Z;
  X.resize(N);
  Y.resize(N);
  Z.resize(N);
  for (int i = 0; i < N; ++i) {
    X[i].resize(N);
    Y[i].resize(N);
    Z[i].resize(N);
  }

  auto z = read_access_host(dynamic_cast<const CahnHilliardState&>(state).c());

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      X[i][j] = i;
      Y[i][j] = j;
      Z[i][j] = z(i, j);
    }
  }

  matplot::figure();
  matplot::contour(X, Y, Z);  // contourf segfaults for some reason?
  matplot::xlabel("X");
  matplot::ylabel("Y");

  matplot::show();
}

std::unique_ptr<Discretization2DCart>
Puppeteer::geomFactory()
{
  const int nhalo = 2;  // should be dependent on equation and order / stencil size
  return std::make_unique<Discretization2DCart>(
      Options::get().domain_resolution(),
      nhalo,
      Options::get().domain_x_begin(),
      Options::get().domain_x_end(),
      Options::get().domain_x_begin());
}


std::unique_ptr<TimeIntegrableRHS>
Puppeteer::rhsFactory(Discretization2DCart& geom)
{
  Logger::get().FatalAssert(
      Options::get().equation_type() == "cahn-hilliard", "Only equation_type=cahn-hilliard currently supported.");

  CahnHilliardParameters params;
  return std::make_unique<CahnHilliard2DFD>(geom, params);
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
Puppeteer::attachSolutionObservers(TimeIntegrator& time_integrator, Discretization2DCart& geom, TimeIntegrableRHS& rhs)
{
  time_integrator.registerObserver(Event::TimeStepComplete, [&]() {
    Logger::get().setSolutionMonitor(
        cahn_hilliard_solution_monitor(rhs, geom, time_integrator.getCurrentSolutionState()));
  });
}
