#pragma once


#include "governing_equations/initial_conditions.hpp"

class VTKSolutionReader;


class InitialConditionsFromFile : public InitialConditions {
 public:
  InitialConditionsFromFile(VTKSolutionReader& reader) : reader_(reader) {}

  virtual ~InitialConditionsFromFile() = default;

  void set(const TimeIntegrableRHS&, SolutionState& state) const override;

 private:
  VTKSolutionReader& reader_;
};
