#pragma once

#include "governing_equations/initial_conditions.hpp"

class CahnHilliardParameters;

class CahnHilliardInitialConditions : public InitialConditions
{
 public:
  CahnHilliardInitialConditions(CahnHilliardParameters& params);
  virtual ~CahnHilliardInitialConditions() = default;

  void set(const TimeIntegrableRHS&, SolutionState& state) const override;

 private:
  const double m_;
  const double min_;
  const double max_;
};
