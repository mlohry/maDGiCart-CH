#pragma once


#include "cahn_hilliard_base.hpp"
#include "cahn_hilliard_parameters.hpp"
#include "cahn_hilliard_state.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"


class CahnHilliard3DFD : public CahnHilliardBase {
 public:
  CahnHilliard3DFD(Discretization3DCart& geom, const CahnHilliardParameters& params);

  void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) override;

  std::unique_ptr<SolutionState> createSolutionState() const override
  {
    return std::make_unique<CahnHilliardState3D>(
        ManagedArrayOwner{}, "state", geom_.ni(), geom_.ni(), geom_.ni(), geom_.nhalo());
  }

  int_t dofsPerEquation() const override { return geom_.nInteriorPoints(); }

  const IndexArray& interiorIndices() const override { return geom_.interiorIndices(); }

  void evalRHSBasic(const SolutionState& flovars, double time, SolutionState& rhs);
  void evalRHSFused(const SolutionState& flovars, double time, SolutionState& rhs);
  void evalRHSFullyFused(const SolutionState& flovars, double time, SolutionState& rhs);

 private:
  Discretization3DCart& geom_;
  const int kernel_variant_;

  std::unique_ptr<ManagedArray3D<real_wp>> laplacian_rhs_term_;
  std::unique_ptr<ManagedArray3D<real_wp>> laplacian_argument_;
  std::unique_ptr<ManagedArray3D<real_wp>> biharmonic_term_;
  std::unique_ptr<ManagedArray3D<real_wp>> linear_term_;


};
