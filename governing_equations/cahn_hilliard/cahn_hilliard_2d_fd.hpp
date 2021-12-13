#pragma once

#include "cahn_hilliard_parameters.hpp"
#include "cahn_hilliard_state.hpp"
#include "data_structures/managed_array_2d.hpp"
#include "governing_equations/time_integrable_rhs.hpp"
#include "logger/logger.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"


class CahnHilliard2DFD : public TimeIntegrableRHS {
 public:
  CahnHilliard2DFD(Discretization2DCart& geom, const CahnHilliardParameters& params);


  void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) override;


  std::unique_ptr<SolutionState> createSolutionState() const override
  {
    return std::make_unique<CahnHilliardState>(ManagedArrayOwner{}, "state", geom_.ni(), geom_.ni(), geom_.nhalo());
  }


  int_t dofsPerEquation() const override { return geom_.nInteriorPoints(); }
  int_t nEquations() const override { return 1; }
  std::vector<std::string> equationNames() const override { return {{"c"}}; }

  void setRandomInitialCondition(SolutionState& state) const;


  std::map<std::string, double> solutionReport(const SolutionState& state, const SolutionState& residual);

  const IndexArray& interiorIndices() const override { return geom_.interiorIndices(); }

  double timeStepLimit(const SolutionState& solution) override { Logger::get().FatalMessage("not implemented"); }

 private:
  Discretization2DCart& geom_;

  const double m_;
  const double sigma_;
  const double eps2_;
  const double initial_value_min_;
  const double initial_value_max_;
};
