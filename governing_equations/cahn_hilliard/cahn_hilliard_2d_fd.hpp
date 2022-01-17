#pragma once

#include "cahn_hilliard_parameters.hpp"
#include "cahn_hilliard_state.hpp"
#include "cahn_hilliard_base.hpp"
#include "data_structures/managed_array_2d.hpp"
#include "governing_equations/time_integrable_rhs.hpp"
#include "spatial_discretization/discretization_2d_cart.hpp"


class CahnHilliard2DFD : public CahnHilliardBase {
 public:
  CahnHilliard2DFD(Discretization2DCart& geom, const CahnHilliardParameters& params);

  void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) override;

  std::unique_ptr<SolutionState> createSolutionState() const override
  {
    return std::make_unique<CahnHilliardState>(ManagedArrayOwner{}, "state", geom_.ni(), geom_.ni(), geom_.nhalo());
  }

  int_t dofsPerEquation() const override { return geom_.nInteriorPoints(); }

  const IndexArray& interiorIndices() const override { return geom_.interiorIndices(); }

 private:
  Discretization2DCart& geom_;

  std::unique_ptr<ManagedArray2D<real_wp>> laplacian_rhs_term_;
  std::unique_ptr<ManagedArray2D<real_wp>> laplacian_argument_;
  std::unique_ptr<ManagedArray2D<real_wp>> biharmonic_term_;
  std::unique_ptr<ManagedArray2D<real_wp>> linear_term_;
};
