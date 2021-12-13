#pragma once

#include <memory>
#include <string>
#include <vector>

#include "program_options/program_options_parser.hpp"

class SolutionState;
class TimeIntegrableRHS;
class Discretization2DCart;
class InitialConditions;
class TimeIntegrator;


class Puppeteer {
 public:
  Puppeteer(const std::vector<std::string>& cmd_line);

  ~Puppeteer();

  void run();

 private:
  std::unique_ptr<Discretization2DCart> geom_;
  std::unique_ptr<TimeIntegrableRHS>    rhs_;
  std::unique_ptr<InitialConditions>    initial_conditions_;
  std::unique_ptr<TimeIntegrator>       time_integrator_;

  static std::unique_ptr<TimeIntegrableRHS>    rhsFactory(Discretization2DCart&);
  static std::unique_ptr<Discretization2DCart> geomFactory();
  static std::unique_ptr<InitialConditions>    initialConditionsFactory();
  static std::unique_ptr<TimeIntegrator>       timeIntegratorFactory(TimeIntegrableRHS&);

  void attachTimeObservers(TimeIntegrator&);
  void attachResidualObservers(TimeIntegrator&, TimeIntegrableRHS&);
  void attachSolutionObservers(TimeIntegrator&, Discretization2DCart&, TimeIntegrableRHS&);

  void plotSolution(const SolutionState&);
};
