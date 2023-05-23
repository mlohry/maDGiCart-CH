#pragma once

#include "governing_equations/initial_conditions.hpp"


class ConstantInitialConditions : public InitialConditions
{
 public:
  ConstantInitialConditions(double value) : value_(value) {}
  virtual ~ConstantInitialConditions() = default;

  void set(const TimeIntegrableRHS&, SolutionState& state) const override;

 private:
  const double value_;
};
