#include "discretization_3d_cart.hpp"

Discretization3DCart::Discretization3DCart(int ni, int nhalo, double xbeg, double xend, double ybeg, double zbeg)
    : ManagedArrayOwner("Discretization3DCart"),
      ni_(ni),
      nhalo_(nhalo),
      ninhalo_(ni + 2 * nhalo),
      dx_((xend - xbeg) / real_wp(ni)),
      interior_indices_(*this, "InteriorIndices", ni * ni * ni),
      periodic_donor_indices_(*this, "PeriodicDonorIndices", (pow(ninhalo_, 3) - pow(ni_, 3))),
      periodic_receiver_indices_(*this, "PeriodicReceiverIndices", (pow(ninhalo_, 3) - pow(ni_, 3))),
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
    auto x = write_access_host(x_coord_);
    auto y = write_access_host(y_coord_);
    auto z = write_access_host(z_coord_);

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
