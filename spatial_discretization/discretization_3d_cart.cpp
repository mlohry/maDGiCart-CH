#include "discretization_3d_cart.hpp"

Discretization3DCart::Discretization3DCart(int ni, int nhalo, double xbeg, double xend, double ybeg, double zbeg)
    : ManagedArrayOwner("Discretization3DCart"),
      ni_(ni),
      nhalo_(nhalo),
      ninhalo_(ni + 2 * nhalo),
      dx_((xend - xbeg) / real_wp(ni)),
      interior_indices_(*this, "InteriorIndices", ni * ni * ni),
      x_coord_(*this, "CoordinateX", ni, ni, ni, nhalo),
      y_coord_(*this, "CoordinateY", ni, ni, ni, nhalo),
      z_coord_(*this, "CoordinateZ", ni, ni, ni, nhalo)
{
  {
    auto idx_list = write_access_host(interior_indices_);
    int  inode    = 0;

    for (int i = 0; i < ni; ++i) {
      for (int j = 0; j < ni; ++j) {
        for (int k = 0; k < ni; ++k) {
          idx_list[inode] = get1DindexFrom3D(i, j, k, nhalo_, ninhalo_, ninhalo_);
          inode++;
        }
      }
    }
  }


  {
    auto          x      = write_access_host(x_coord_);
    auto          y      = write_access_host(y_coord_);
    auto          z      = write_access_host(z_coord_);
    const real_wp xc_beg = xbeg + 0.5 * dx_;
    const real_wp yc_beg = ybeg + 0.5 * dx_;
    const real_wp zc_beg = zbeg + 0.5 * dx_;
    for (int i = -nhalo; i < ni + nhalo; ++i) {
      for (int j = -nhalo; j < ni + nhalo; ++j) {
        for (int k = -nhalo; k < ni + nhalo; ++k) {
          x(i, j, k) = xc_beg + real_wp(i) * dx_;
          y(i, j, k) = yc_beg + real_wp(j) * dx_;
          z(i, j, k) = zc_beg + real_wp(k) * dx_;
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
Discretization3DCart::applyPeriodicBoundaryConditions(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto      state  = read_write_access(state_in);
  const int nhalos = nhalo();
  const int nis    = ni();
  const int njs    = ni();
  const int nks    = ni();

  // x/i halos
  maDGForAll(j, 0, njs, {
    for (int k = 0; k < nks; ++k) {
      for (int i = -nhalos; i < 0; ++i) {
        state(i, j, k) = state(i + nis, j, k);
      }
      for (int i = nis; i < nis + nhalos; ++i) {
        state(i, j, k) = state(i - nis, j, k);
      }
    }
  });


  // y/j halos
  maDGForAll(i, 0, nis, {
    for (int k = 0; k < nks; ++k) {
      for (int j = -nhalos; j < 0; ++j) {
        state(i, j, k) = state(i, j + njs, k);
      }
      for (int j = njs; j < njs + nhalos; ++j) {
        state(i, j, k) = state(i, j - njs, k);
      }
    }
  });


  // z/k halos
  maDGForAll(i, 0, nis, {
    for (int j = 0; j < njs; ++j) {
      for (int k = -nhalos; k < 0; ++k) {
        state(i, j, k) = state(i, j, k + nks);
      }
      for (int k = nks; k < nks + nhalos; ++k) {
        state(i, j, k) = state(i, j, k - nks);
      }
    }
  });

//  // corner diagonals aligned with k
//  maDGForAll(k, 0, nks, {
//    for (int i = -nhalos; i < 0; ++i) {
//      for (int j = -nhalos; j < 0; ++j) {
//        state(i, j, k) = state(i + nis, j + njs, k);
//      }
//      for (int j = njs; j < njs + nhalos; ++j) {
//        state(i, j, k) = state(i + nis, j - njs, k);
//      }
//    }
//    for (int i = nis; i < nis + nhalos; ++i) {
//      for (int j = -nhalos; j < 0; ++j) {
//        state(i, j, k) = state(i - nis, j + njs, k);
//      }
//      for (int j = njs; j < njs + nhalos; ++j) {
//        state(i, j, k) = state(i - nis, j - nis, k);
//      }
//    }
//  });
}
