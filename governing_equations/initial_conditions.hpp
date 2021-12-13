#pragma once

#include "typedefs.hpp"

class TimeIntegrableRHS;
class SolutionState;

class InitialConditions
{
 public:
  virtual ~InitialConditions() = default;

  virtual void set(const TimeIntegrableRHS&, SolutionState& state) const = 0;
};
