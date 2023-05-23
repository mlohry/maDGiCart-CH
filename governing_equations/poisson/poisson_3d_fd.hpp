#pragma once


#include "data_structures/scalar_solution_state.hpp"
#include "governing_equations/time_integrable_rhs.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"


class Poisson3DFD : public TimeIntegrableRHS {
 public:
  Poisson3DFD(Discretization3DCart& geom);

  int_t nEquations() const override { return 1; }

  std::vector<std::string> equationNames() const override { return {{"c"}}; }

  void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) override;

  void addSourceTerm(ScalarSolutionState3D& rhs);

  std::unique_ptr<SolutionState> createSolutionState() const override
  {
    return std::make_unique<ScalarSolutionState3D>(
        ManagedArrayOwner{}, "state", geom_.ni(), geom_.nj(), geom_.nk(), geom_.nhalo());
  }

  const SpatialDiscretization* hasGeometry() const override { return &geom_; }

  int_t dofsPerEquation() const override { return geom_.nInteriorPoints(); }

  const IndexArray& interiorIndices() const override { return geom_.interiorIndices(); }

  void applyDirichletBoundaryConditions(ScalarSolutionState3D& state);

  std::unique_ptr<TimeIntegrableRHS> clone(SpatialDiscretization& geom) const override;

  double cflScale() const override
  {
    const double dx = (geom_.domain().xend - geom_.domain().xbeg) / double(geom_.domain().nx);
    return dx * dx;
  }

 private:
  Discretization3DCart& geom_;
};
