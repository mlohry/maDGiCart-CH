#pragma once

#include "data_structures/managed_array_3d.hpp"
#include "spatial_discretization.hpp"

class Discretization3DCart : public SpatialDiscretization {
 public:

  Discretization3DCart(const CartesianDomainDefinition& domain);

  ~Discretization3DCart() = default;

  const IndexArray& interiorIndices() const { return interior_indices_; }
  //  const IndexArray& boundaryIndices() const { return boundary_indices_; }

  std::unique_ptr<ManagedArray3D<real_wp>> createRealArray() const;

  int    ni() const { return ni_; }
  int    nhalo() const { return nhalo_; }
  double dx() const { return dx_; }
  int    nInteriorPoints() const { return interior_indices_.size(); }

  void applyPeriodicBoundaryConditions(ManagedArray3D<real_wp>& state) const;
  void applyNeumannBoundaryConditions(ManagedArray3D<real_wp>& state) const;


  void laplacian(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del2state_out) const;
  void biharmonic(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del4state_out) const;

  const ManagedArray3DOwning<real_wp>& x() const { return x_coord_; }
  const ManagedArray3DOwning<real_wp>& y() const { return y_coord_; }
  const ManagedArray3DOwning<real_wp>& z() const { return z_coord_; }

  const ManagedArray3DOwning<real_wp>& xvertex() const { return x_vertex_coord_; }
  const ManagedArray3DOwning<real_wp>& yvertex() const { return y_vertex_coord_; }
  const ManagedArray3DOwning<real_wp>& zvertex() const { return z_vertex_coord_; }

  void applyNeumannBCX(ManagedArray3D<real_wp>& state) const;
  void applyNeumannBCY(ManagedArray3D<real_wp>& state) const;
  void applyNeumannBCZ(ManagedArray3D<real_wp>& state) const;

 private:
  const int    ni_;
  const int    nhalo_;
  const int    ninhalo_;
  const double dx_;

  IndexArray interior_indices_;

  IndexArray periodic_donor_indices_;
  IndexArray periodic_receiver_indices_;

  IndexArray xmin_indices_;
  IndexArray xmax_indices_;
  IndexArray ymin_indices_;
  IndexArray ymax_indices_;
  IndexArray zmin_indices_;
  IndexArray zmax_indices_;

  ManagedArray3DOwning<real_wp> x_coord_;
  ManagedArray3DOwning<real_wp> y_coord_;
  ManagedArray3DOwning<real_wp> z_coord_;

  ManagedArray3DOwning<real_wp> x_vertex_coord_;
  ManagedArray3DOwning<real_wp> y_vertex_coord_;
  ManagedArray3DOwning<real_wp> z_vertex_coord_;
};
