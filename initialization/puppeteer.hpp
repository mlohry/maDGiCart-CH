#pragma once

#include <memory>
#include <string>
#include <vector>

#include "program_options/program_options_parser.hpp"

class SolutionState;
class TimeIntegrableRHS;
class SpatialDiscretization;
class InitialConditions;
class TimeIntegrator;
class VTKSolutionReader;


class Puppeteer {
 public:
  Puppeteer(const std::vector<std::string>& cmd_line);

  ~Puppeteer();

  void run();

 private:
  std::unique_ptr<SpatialDiscretization> geom_;
  std::unique_ptr<TimeIntegrableRHS>     rhs_;
  std::unique_ptr<InitialConditions>     initial_conditions_;
  std::unique_ptr<TimeIntegrator>        time_integrator_;
  std::unique_ptr<VTKSolutionReader>     solution_reader_;

  std::string file_output_string_;

  static std::unique_ptr<TimeIntegrableRHS> rhsFactory(SpatialDiscretization&);
  std::unique_ptr<SpatialDiscretization>    geomFactory();
  std::unique_ptr<InitialConditions>        initialConditionsFactory();
  static std::unique_ptr<TimeIntegrator>    timeIntegratorFactory(TimeIntegrableRHS&);

  void attachTimeObservers(TimeIntegrator&);
  void attachResidualObservers(TimeIntegrator&, TimeIntegrableRHS&);
  void attachSolutionObservers(TimeIntegrator&, SpatialDiscretization&, TimeIntegrableRHS&);
};
