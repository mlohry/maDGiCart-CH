#include "puppeteer.hpp"

#include "file_io/vtk_solution_reader.hpp"
#include "file_io/vtk_solution_writer.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_2d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_3d_fd.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_initial_conditions.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_parameters.hpp"
#include "governing_equations/cahn_hilliard/cahn_hilliard_solution_monitor.hpp"
#include "governing_equations/poisson/poisson_3d_fd.hpp"
#include "governing_equations/initial_conditions.hpp"
#include "governing_equations/constant_initial_conditions.hpp"
#include "governing_equations/initial_conditions_from_file.hpp"
#include "governing_equations/residual_calculator.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"
#include "time_stepping/time_integrator.hpp"
#include "time_stepping/time_integrator_options.hpp"
#include "utils/registry.hpp"

#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <iomanip>


Puppeteer::Puppeteer(const std::vector<std::string>& cmd_line)
{
  ProgramOptionsParser parser;
  parser.parseInputOptions(cmd_line);
  programOptionsSanityCheck();

  Logger::get().initLogFile(Options::get().log_file());
  Logger::get().InfoMessage(parser.getConfigurationFileTemplate());

  const auto init_cond_file = Options::get().initial_condition_file();
  if (init_cond_file != "none") {
    solution_reader_ = std::make_unique<VTKSolutionReader>(init_cond_file);
    Logger::get().InfoMessage("Using initial conditions file " + init_cond_file);
  }

  file_output_string_ = [&]() {
    if (Options::get().solution_output_file().empty()) {
      return Options::get().log_file();
    }
    return Options::get().solution_output_file();
  }();

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

  time_integrator_->setNotifyFrequency(Event::TimeStepComplete, Options::get().log_frequency());
  time_integrator_->setNotifyFrequency(Event::SolutionUpdate, Options::get().save_every());
}


Puppeteer::~Puppeteer() = default;


void
Puppeteer::run()
{
  profile();
  auto run_timer = Logger::get().timer("Starting solver");

  time_integrator_->solve(*rhs_, *initial_conditions_);

  write_solution_to_vtk(
      time_integrator_->getIterationStatus(),
      time_integrator_->getCurrentSolutionState(),
      *geom_,
      file_output_string_ + "_final");

  Logger::get().InfoMessage(Profiler::get().finalize());
}

namespace {
BCType
parseBC(const std::string bcstr)
{
  if (bcstr == "periodic") {
    return BCType::Periodic;
  }
  if (bcstr == "neumann") {
    return BCType::Neumann;
  }
  if (bcstr == "dirichlet") {
    return BCType::Dirichlet;
  }
  Logger::get().FatalMessage("BC type string " + bcstr + " not recognized.");
  abort();
}
}  // namespace


std::unique_ptr<SpatialDiscretization>
Puppeteer::geomFactory()
{
  CartesianDomainDefinition domain;

  if (solution_reader_) {
    domain = solution_reader_->getCartesianDomain();
  }
  else {
    domain.nx   = Options::get().domain_resolution_x();
    domain.ny   = Options::get().domain_resolution_y();
    domain.nz   = Options::get().domain_resolution_z();
    domain.xbeg = Options::get().domain_x_begin();
    domain.xend = Options::get().domain_x_end();
    domain.ybeg = Options::get().domain_y_begin();
    domain.yend = Options::get().domain_y_end();
    domain.zbeg = Options::get().domain_z_begin();
    domain.zend = Options::get().domain_z_end();
  }

  domain.nhalo = 2;  // should be dependent on equation and order / stencil size

  domain.xbc = parseBC(Options::get().bc_x());
  domain.ybc = parseBC(Options::get().bc_y());
  domain.zbc = parseBC(Options::get().bc_z());

  const int dim = Options::get().dimension();
  switch (dim) {
    case 2:
      return std::make_unique<Discretization2DCart>(domain);
    case 3:
      return std::make_unique<Discretization3DCart>(domain);
    default:
      Logger::get().FatalMessage("Only dimensions 2 and 3 supported.");
  }
  return nullptr;
}


std::unique_ptr<TimeIntegrableRHS>
Puppeteer::rhsFactory(SpatialDiscretization& geom)
{
  if (Options::get().equation_type() == "cahn-hilliard") {

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
  }
  else if (Options::get().equation_type() == "poisson") {
    return std::make_unique<Poisson3DFD>(dynamic_cast<Discretization3DCart&>(geom));
  }
  else {
    Logger::get().FatalMessage("Only cahn-hilliard and poission equation_type supported.");
  }



  return nullptr;
}


std::unique_ptr<InitialConditions>
Puppeteer::initialConditionsFactory()
{
  if (solution_reader_) {
    return std::make_unique<InitialConditionsFromFile>(*solution_reader_);
  }
  else {
    if (Options::get().equation_type() == "cahn-hilliard"){
      CahnHilliardParameters params;
      return std::make_unique<CahnHilliardInitialConditions>(params);
    }
    else if (Options::get().equation_type() == "poisson"){
      return std::make_unique<ConstantInitialConditions>(0.0);
    }
    else {
      Logger::get().FatalMessage("Only cahn-hilliard and poisson supported for initialConditionsFactory");
    }
  }
  return nullptr;
}


std::unique_ptr<TimeIntegrator>
Puppeteer::timeIntegratorFactory(TimeIntegrableRHS& rhs)
{
  TimeIntegratorOptions time_opts = getTimeOptionsUsingGlobalOptions();
  if (solution_reader_) {
    solution_reader_->setInitialTimeStep(time_opts);
  }
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
      break;
    case 3:
      time_integrator.registerObserver(Event::TimeStepComplete, [&]() {
        Logger::get().setSolutionMonitor(cahn_hilliard_solution_monitor(
            rhs, dynamic_cast<Discretization3DCart&>(geom), time_integrator.getCurrentSolutionState()));
      });
      break;
    default:
      Logger::get().FatalMessage("Only dimensions 2 and 3 supported.");
  }

  const int save_every = Options::get().save_every();

  if (save_every) {
    time_integrator.registerObserver(Event::SolutionUpdate, [&]() {
      write_solution_to_vtk(
          time_integrator_->getIterationStatus(),
          time_integrator_->getCurrentSolutionState(),
          *geom_,
          file_output_string_);
    });
  }
}


void
Puppeteer::programOptionsSanityCheck() const
{
  //  if (Options::get().dimension() == 3) {
  //    if (Options::get().bc_x() != "periodic" || Options::get().bc_y() != "periodic" ||
  //        Options::get().bc_z() != "periodic") {
  //      if (Options::get().time_integrator().compare(0, 5, "petsc")){
  //        Logger::get().FatalMessage("petsc implicit solvers only support triply-periodic case");
  //      }
  //      if (Options::get().kernel_variant() != 0){
  //        Logger::get().FatalMessage("kernel_variant only supported for triply-periodic case");
  //      }
  //    }
  //  }
}
