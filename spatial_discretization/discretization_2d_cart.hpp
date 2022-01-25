#pragma once

#include "data_structures/managed_array_2d.hpp"
#include "spatial_discretization.hpp"


class Discretization2DCart : public SpatialDiscretization {
 public:
  Discretization2DCart(int ni, int nhalo, double xbeg, double xend, double ybeg);

  Discretization2DCart(const CartesianDomainDefinition& domain, int nhalo)
      : Discretization2DCart(domain.nx, nhalo, domain.xbeg, domain.xend, domain.ybeg)
  {
  }


  ~Discretization2DCart() = default;

  const IndexArray& interiorIndices() const { return interior_indices_; }
  const IndexArray& boundaryIndices() const { return boundary_indices_; }

  std::unique_ptr<ManagedArray2D<real_wp>> createRealArray() const;

  int    ni() const { return ni_; }
  int    nhalo() const { return nhalo_; }
  double dx() const { return dx_; }
  int    nInteriorPoints() const { return ni() * ni(); }

  void applyPeriodicBoundaryConditions(ManagedArray2D<real_wp>& state) const;

  void printState(const ManagedArray2D<real_wp>& state) const;

  void laplacian(const ManagedArray2D<real_wp>& state_in, ManagedArray2D<real_wp>& del2state_out) const;
  void biharmonic(const ManagedArray2D<real_wp>& state_in, ManagedArray2D<real_wp>& del4state_out) const;

  const ManagedArray2DOwning<real_wp>& x() const { return x_coord_; }
  const ManagedArray2DOwning<real_wp>& y() const { return y_coord_; }

  const ManagedArray2DOwning<real_wp>& xvertex() const { return x_vertex_coord_; }
  const ManagedArray2DOwning<real_wp>& yvertex() const { return y_vertex_coord_; }

 private:
  const int    ni_;
  const int    nhalo_;
  const int    ninhalo_;
  const double dx_;

  IndexArray                    interior_indices_;
  IndexArray                    boundary_indices_;
  ManagedArray2DOwning<real_wp> x_coord_;
  ManagedArray2DOwning<real_wp> y_coord_;

  ManagedArray2DOwning<real_wp> x_vertex_coord_;
  ManagedArray2DOwning<real_wp> y_vertex_coord_;
};
