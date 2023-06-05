#pragma once


#include "cahn_hilliard_base.hpp"
#include "cahn_hilliard_parameters.hpp"
#include "data_structures/scalar_solution_state.hpp"
#include "spatial_discretization/discretization_3d_cart.hpp"


class CahnHilliard3DFD : public CahnHilliardBase {
 public:
  CahnHilliard3DFD(Discretization3DCart& geom, const CahnHilliardParameters& params);

  void evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs) override;

  void evalJacobian(const SolutionState& flovars, double time, CSRMatrix& J) override;

  std::unique_ptr<SolutionState> createSolutionState() const override
  {
    return std::make_unique<ScalarSolutionState3D>(
        ManagedArrayOwner{}, "state", geom_.ni(), geom_.nj(), geom_.nk(), geom_.nhalo());
  }

  int_t dofsPerEquation() const override { return geom_.nInteriorPoints(); }

  const IndexArray& interiorIndices() const override { return geom_.interiorIndices(); }

  std::unique_ptr<CSRMatrix> createSparseMatrix() const override;
  std::vector<int>           nNonZerosPerRow() const override;

  void evalRHSBasic(const SolutionState& flovars, double time, SolutionState& rhs);
  void evalRHSFused(const SolutionState& flovars, double time, SolutionState& rhs);
  void evalRHSFullyFused(const SolutionState& flovars, double time, SolutionState& rhs);

  const SpatialDiscretization* hasGeometry() const override { return &geom_; }

  std::unique_ptr<TimeIntegrableRHS> clone(SpatialDiscretization&) const override;

  double cflScale() const override
  {
    const double dx = (geom_.domain().xend - geom_.domain().xbeg) / double(geom_.domain().nx);
    return dx * dx * dx * dx;
  }


 private:
  Discretization3DCart& geom_;
  const int             kernel_variant_;

  std::unique_ptr<ManagedArray3D<real_wp>> laplacian_rhs_term_;
  std::unique_ptr<ManagedArray3D<real_wp>> laplacian_argument_;
  std::unique_ptr<ManagedArray3D<real_wp>> biharmonic_term_;
  std::unique_ptr<ManagedArray3D<real_wp>> linear_term_;
};
