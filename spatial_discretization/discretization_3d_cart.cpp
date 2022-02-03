#include "discretization_3d_cart.hpp"



Discretization3DCart::Discretization3DCart(const CartesianDomainDefinition& domain)
    : SpatialDiscretization("Discretization3DCart", domain),
      ni_(domain.nx),
      nhalo_(domain.nhalo),
      ninhalo_(ni_ + 2 * nhalo_),
      dx_((domain.xend - domain.xbeg) / real_wp(ni_)),
      interior_indices_(*this, "InteriorIndices", ni_ * ni_ * ni_),
      periodic_donor_indices_(*this, "PeriodicDonorIndices", (pow(ninhalo_, 3) - pow(ni_, 3))),
      periodic_receiver_indices_(*this, "PeriodicReceiverIndices", (pow(ninhalo_, 3) - pow(ni_, 3))),
      xmin_indices_(*this, "XMinIndices", ninhalo_ * ninhalo_ * nhalo_),
      xmax_indices_(*this, "XMaxIndices", ninhalo_ * ninhalo_ * nhalo_),
      ymin_indices_(*this, "YMinIndices", ninhalo_ * ninhalo_ * nhalo_),
      ymax_indices_(*this, "YMaxIndices", ninhalo_ * ninhalo_ * nhalo_),
      zmin_indices_(*this, "ZMinIndices", ninhalo_ * ninhalo_ * nhalo_),
      zmax_indices_(*this, "ZMaxIndices", ninhalo_ * ninhalo_ * nhalo_),
      x_coord_(*this, "CoordinateX", ni_, ni_, ni_, nhalo_),
      y_coord_(*this, "CoordinateY", ni_, ni_, ni_, nhalo_),
      z_coord_(*this, "CoordinateZ", ni_, ni_, ni_, nhalo_),
      x_vertex_coord_(*this, "VertexCoordinateX", ni_ + 1, ni_ + 1, ni_ + 1, nhalo_),
      y_vertex_coord_(*this, "VertexCoordinateY", ni_ + 1, ni_ + 1, ni_ + 1, nhalo_),
      z_vertex_coord_(*this, "VertexCoordinateZ", ni_ + 1, ni_ + 1, ni_ + 1, nhalo_)
{
  {
    auto idx_list = write_access_host(interior_indices_);
    int  inode    = 0;

    for (int i = 0; i < ni_; ++i) {
      for (int j = 0; j < ni_; ++j) {
        for (int k = 0; k < ni_; ++k) {
          idx_list[inode] = get1DindexFrom3D(i, j, k, nhalo_, ninhalo_, ninhalo_);
          inode++;
        }
      }
    }
  }

  {
    auto xv = write_access_host(x_vertex_coord_);
    auto yv = write_access_host(y_vertex_coord_);
    auto zv = write_access_host(z_vertex_coord_);

    for (int i = -nhalo_; i < ni_ + nhalo_ + 1; ++i) {
      for (int j = -nhalo_; j < ni_ + nhalo_ + 1; ++j) {
        for (int k = -nhalo_; k < ni_ + nhalo_ + 1; ++k) {
          xv(i, j, k) = domain.xbeg + dx_ * real_wp(i);
          yv(i, j, k) = domain.ybeg + dx_ * real_wp(j);
          zv(i, j, k) = domain.zbeg + dx_ * real_wp(k);
        }
      }
    }
  }


  {
    auto xv = read_access_host(x_vertex_coord_);
    auto yv = read_access_host(y_vertex_coord_);
    auto zv = read_access_host(z_vertex_coord_);
    auto x  = write_access_host(x_coord_);
    auto y  = write_access_host(y_coord_);
    auto z  = write_access_host(z_coord_);

    for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
      for (int j = -nhalo_; j < ni_ + nhalo_; ++j) {
        for (int k = -nhalo_; k < ni_ + nhalo_; ++k) {
          x(i, j, k) = 0.5 * (xv(i, j, k) + xv(i + 1, j, k));
          y(i, j, k) = 0.5 * (yv(i, j, k) + yv(i, j + 1, k));
          z(i, j, k) = 0.5 * (zv(i, j, k) + zv(i, j, k + 1));
        }
      }
    }
  }

  {
    int  inode         = 0;
    auto donor_list    = write_access_host(periodic_donor_indices_);
    auto receiver_list = write_access_host(periodic_receiver_indices_);

    for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
      for (int j = -nhalo_; j < ni_ + nhalo_; ++j) {
        for (int k = -nhalo_; k < ni_ + nhalo_; ++k) {
          bool is_halo_cell = [&]() {
            if (i < 0 || i >= ni_ || j < 0 || j >= ni_ || k < 0 || k >= ni_) {
              return true;
            }
            return false;
          }();

          if (is_halo_cell) {
            receiver_list[inode] = get1DindexFrom3D(i, j, k, nhalo_, ninhalo_, ninhalo_);

            const int idonor = [&]() {
              if (i < 0) {
                return i + ni_;
              }
              if (i >= ni_) {
                return i - ni_;
              }
              return i;
            }();

            const int jdonor = [&]() {
              if (j < 0) {
                return j + ni_;
              }
              if (j >= ni_) {
                return j - ni_;
              }
              return j;
            }();

            const int kdonor = [&]() {
              if (k < 0) {
                return k + ni_;
              }
              if (k >= ni_) {
                return k - ni_;
              }
              return k;
            }();

            donor_list[inode] = get1DindexFrom3D(idonor, jdonor, kdonor, nhalo_, ninhalo_, ninhalo_);

            inode++;
          }
        }
      }
    }
  }
}


std::unique_ptr<ManagedArray3D<real_wp>>
Discretization3DCart::createRealArray() const
{
  return std::make_unique<ManagedArray3DOwning<real_wp>>(*this, "temp", ni_, ni_, ni_, nhalo_);
}


void
Discretization3DCart::laplacian(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del2state_out) const
{
  auto state     = read_access(state_in);
  auto del2state = write_access(del2state_out);
  auto idx       = read_access(interiorIndices());

  const real_wp dx2 = dx() * dx();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    int k;
    state.getIJK(idx[ii], i, j, k);

    del2state(i, j, k) = (state(i - 1, j, k) + state(i + 1, j, k) + state(i, j - 1, k) + state(i, j + 1, k) +
                          state(i, j, k - 1) + state(i, j, k + 1) - real_wp(6) * state(i, j, k)) /
                         dx2;
  });
}


void
Discretization3DCart::biharmonic(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del4state_out) const
{
  ///
  ///  High accuracy solution of three-dimensional biharmonic equations, Altas et al
  /// https://www.irisa.fr/sage/jocelyne/publis/2000/num-alg-2002.pdf
  ///
  /// https://arxiv.org/pdf/1901.05118.pdf
  assert(nhalo_ >= 2);

  auto f     = read_access(state_in);
  auto del4f = write_access(del4state_out);
  auto idx   = read_access(interiorIndices());

  const real_wp dx4 = dx() * dx() * dx() * dx();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    int k;
    f.getIJK(idx[ii], i, j, k);

    del4f(i, j, k) =
        (real_wp(42) * f(i, j, k) -
         real_wp(12) *
             (f(i + 1, j, k) + f(i, j + 1, k) + f(i, j, k + 1) + f(i - 1, j, k) + f(i, j - 1, k) + f(i, j, k - 1)) +
         f(i + 2, j, k) + f(i, j + 2, k) + f(i, j, k + 2) + f(i - 2, j, k) + f(i, j - 2, k) + f(i, j, k - 2) +
         real_wp(2) * (f(i - 1, j - 1, k) + f(i - 1, j + 1, k) + f(i + 1, j - 1, k) + f(i + 1, j + 1, k) +
                       f(i - 1, j, k - 1) + f(i + 1, j, k - 1) + f(i, j - 1, k - 1) + f(i, j + 1, k - 1) +
                       f(i - 1, j, k + 1) + f(i + 1, j, k + 1) + f(i, j - 1, k + 1) + f(i, j + 1, k + 1))) /
        dx4;
  });
}


void
Discretization3DCart::applyPeriodicBoundaryConditions(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state    = read_write_access(state_in.asArray());
  auto donor    = read_access(periodic_donor_indices_);
  auto receiver = read_access(periodic_receiver_indices_);

  maDGForAll(i, 0, donor.size(), { state[receiver[i]] = state[donor[i]]; });
}


void
Discretization3DCart::applyNeumannBoundaryConditions(ManagedArray3D<real_wp>& state) const
{
  if (domain().xbc == BCType::Neumann){
    applyNeumannBCX(state);
  }
  if (domain().ybc == BCType::Neumann){
    applyNeumannBCY(state);
  }
  if (domain().zbc == BCType::Neumann){
    applyNeumannBCZ(state);
  }
}

void
Discretization3DCart::applyNeumannBCX(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto ninhalo = ninhalo_;
  const auto nhalo = nhalo_;
  const auto ni = ni_;

  maDGForAll(j, -nhalo, ninhalo, {
    for (int k = -nhalo; k < ninhalo; ++k) {
      for (int i = -nhalo; i < 0; ++i){
        state(i, j, k) = state(0, j, k);
      }
      for (int i = ni; i < ninhalo; ++i){
        state(i, j, k) = state(ni-1, j, k);
      }
    }
  });
}

void
Discretization3DCart::applyNeumannBCY(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto ninhalo = ninhalo_;
  const auto nhalo = nhalo_;
  const auto ni = ni_;

  maDGForAll(i, -nhalo, ninhalo, {
    for (int k = -nhalo; k < ninhalo; ++k) {
      for (int j = -nhalo; j < 0; ++j){
        state(i, j, k) = state(i, 0, k);
      }
      for (int j = ni; j < ninhalo; ++j){
        state(i, j, k) = state(i, ni-1, k);
      }
    }
  });
}

void
Discretization3DCart::applyNeumannBCZ(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto ninhalo = ninhalo_;
  const auto nhalo = nhalo_;
  const auto ni = ni_;

  maDGForAll(i, -nhalo, ninhalo, {
    for (int j = -nhalo; j < ninhalo; ++j) {
      for (int k = -nhalo; k < 0; ++k){
        state(i, j, k) = state(i, j, 0);
      }
      for (int k = ni; k < ninhalo; ++k){
        state(i, j, k) = state(i, j, ni-1);
      }
    }
  });
}
