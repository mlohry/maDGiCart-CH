#pragma once

#include "data_structures/managed_array_3d.hpp"
#include "spatial_discretization.hpp"

class Discretization3DCart : public SpatialDiscretization {
 public:
  Discretization3DCart(int ni, int nhalo, double xbeg, double xend, double ybeg, double zbeg);

  ~Discretization3DCart() = default;

  const IndexArray& interiorIndices() const { return interior_indices_; }
//  const IndexArray& boundaryIndices() const { return boundary_indices_; }

  std::unique_ptr<ManagedArray3D<real_wp>> createRealArray() const;

  int ni() const { return ni_; }
  int nhalo() const { return nhalo_; }
  double dx() const { return dx_; }
  int nInteriorPoints() const { return interior_indices_.size(); }

  void applyPeriodicBoundaryConditions(ManagedArray3D<real_wp>& state) const;


  void laplacian(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del2state_out) const;
  void biharmonic(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del4state_out) const;

  const ManagedArray3DOwning<real_wp>& x() const { return x_coord_; }
  const ManagedArray3DOwning<real_wp>& y() const { return y_coord_; }
  const ManagedArray3DOwning<real_wp>& z() const { return z_coord_; }

 private:
  const int ni_;
  const int nhalo_;
  const int ninhalo_;
  const double dx_;

  IndexArray interior_indices_;

  IndexArray periodic_donor_indices_;
  IndexArray periodic_receiver_indices_;

  ManagedArray3DOwning<real_wp> x_coord_;
  ManagedArray3DOwning<real_wp> y_coord_;
  ManagedArray3DOwning<real_wp> z_coord_;
};
