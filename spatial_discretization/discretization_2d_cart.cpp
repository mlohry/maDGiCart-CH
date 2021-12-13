#include "discretization_2d_cart.hpp"

#include <memory>

Discretization2DCart::Discretization2DCart(
    int    ni,
    int    nhalo,
    double xbeg,
    double xend,
    double ybeg)
    : ManagedArrayOwner("Discretization2DCart"),
      ni_(ni),
      nhalo_(nhalo),
      ninhalo_(ni + 2 * nhalo),
      dx_((xend - xbeg) / real_wp(ni)),
      interior_indices_(*this, "InteriorIndices", ni * ni),
      boundary_indices_(*this, "BoundaryIndices", (ni + 2 * nhalo) * nhalo_ * 2 + ni * nhalo_ * 2),
      x_coord_(*this, "CoordinateX", ni, ni, nhalo),
      y_coord_(*this, "CoordinateY", ni, ni, nhalo)
{
  {
    auto idx_list = write_access_host(interior_indices_);
    int  inode    = 0;

    for (int i = 0; i < ni; ++i) {
      for (int j = 0; j < ni; ++j) {
        idx_list[inode] = get1Dindex(i, j, nhalo_, ninhalo_);
        inode++;
      }
    }
  }

  {
    int ibcnode = 0;
    // bottom rows including corners
    auto bc_list = write_access_host(boundary_indices_);
    for (int j = -nhalo; j < 0; ++j) {
      for (int i = -nhalo; i < ni + nhalo; ++i) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // top rows including corners
    for (int j = ni; j < ni + nhalo; ++j) {
      for (int i = -nhalo; i < ni + nhalo; ++i) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // left excluding corners
    for (int i = -nhalo; i < 0; ++i) {
      for (int j = 0; j < ni; ++j) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // right excluding corners
    for (int i = ni; i < ni + nhalo; ++i) {
      for (int j = 0; j < ni; ++j) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }
  }

  auto x = write_access_host(x_coord_);
  auto y = write_access_host(y_coord_);
  const real_wp xc_beg = xbeg+0.5*dx_;
  const real_wp yc_beg = ybeg+0.5*dx_;
//  std::cout << "dy: " << dy_ << " " << yc_beg << "\n";
//  std::cout << "dx: " << dx_ << " " << xc_beg << "\n";

  for (int i = -nhalo; i < ni+nhalo; ++i){
    for (int j = -nhalo; j < ni+nhalo; ++j){
      x(i,j) = xc_beg + real_wp(i)*dx_;
      y(i,j) = yc_beg + real_wp(j)*dx_;
//      std::cout << x(i,j) << " " << y(i,j) << "\n";
    }
  }
}


void
Discretization2DCart::printState(const ManagedArray2D<real_wp>& state_in) const
{
  auto state = read_access_host(state_in);
//  std::cout << "njhalo_-1 " << njhalo_ - 1 << " -nhalo" << -nhalo_ << "\n";
//  std::cout << "-nhalo_ " << -nhalo_ << " nihalo_ " << nihalo_ << "\n";
  for (int j = ni_ + nhalo_ - 1; j >= -nhalo_; --j) {
    for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
      std::cout << std::scientific << state(i, j) << " ";
    }
    std::cout << std::endl;
  }


  // auto state = read_access_host(state_in.asArray());
  //  maDGForAllHost(i, 0, state.size(), {
  //    std::cout << state[i] << "\n";
  //  })
  //
  //  for (int j = njhalo_-1; j >= -nhalo_; --j){
  //    for (int i = -nhalo_; i < nihalo_; ++i){
  //      std::cout << std::scientific << state(i,j) << " ";
  //    }
  //    std::cout << std::endl;
  //  }
}

void
Discretization2DCart::applyPeriodicBoundaryConditions(ManagedArray2D<real_wp>& state_in) const
{
  profile();
  auto      state  = read_write_access(state_in);
  const int nhalos = nhalo();
  const int jstart = nhalos;
  const int jend   = nhalos + ni();
  const int istart = nhalos;
  const int nis    = ni();
  const int njs    = ni();

  // x halos
  maDGForAll(j, 0, njs, {
    for (int i = -nhalos; i < 0; ++i) {
      state(i, j) = state(i + nis, j);
    }
    for (int i = nis; i < nis + nhalos; ++i) {
      state(i, j) = state(i - nis, j);
    }
  });

  // y halos
  maDGForAll(i, 0, nis, {
    for (int j = -nhalos; j < 0; ++j) {
      state(i, j) = state(i, j + njs);
    }
    for (int j = njs; j < njs + nhalos; ++j) {
      state(i, j) = state(i, j - njs);
    }
  });


  maDGForAll(i, -nhalos, 0, {
    for (int j = -nhalos; j < 0; ++j) {  // sw
      state(i, j) = state(i + nis, j + njs);
    }
    for (int j = njs; j < njs + nhalos; ++j) {  // nw
      state(i, j) = state(i + nis, j - njs);
    }
  });

  maDGForAll(i, nis, nis + nhalos, {
    for (int j = -nhalos; j < 0; ++j) {  // se
      state(i, j) = state(i - nis, j + njs);
    }
    for (int j = njs; j < njs + nhalos; ++j) {  // ne
      state(i, j) = state(i - nis, j - njs);
    }
  });
}


void
Discretization2DCart::laplacian(const ManagedArray2D<real_wp>& state_in, ManagedArray2D<real_wp>& del2state_out) const
{
  auto state     = read_access(state_in);
  auto del2state = write_access(del2state_out);
  auto idx       = read_access(interiorIndices());

  const real_wp dx2 = dx() * dx();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    state.getIJ(idx[ii], i, j);

    del2state(i, j) = (state(i - 1, j) - real_wp(2) * state(i, j) + state(i + 1, j)) / dx2 +
                      (state(i, j - 1) - real_wp(2) * state(i, j) + state(i, j + 1)) / dx2;
  });
}


void
Discretization2DCart::biharmonic(const ManagedArray2D<real_wp>& state_in, ManagedArray2D<real_wp>& del4state_out) const
{
  assert(nhalo_ >= 2);

  auto f     = read_access(state_in);
  auto del4f = write_access(del4state_out);
  auto idx   = read_access(interiorIndices());

  const real_wp dx4 = dx() * dx() * dx() * dx();

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    f.getIJ(idx[ii], i, j);

    del4f(i, j) = (f(i, j + 2) + 2.0 * f(i - 1, j + 1) - 8.0 * f(i, j + 1) + 2.0 * f(i + 1, j + 1) + f(i - 2, j) -
                  8.0 * f(i - 1, j) + 20.0 * f(i, j) - 8.0 * f(i + 1, j) + f(i + 2, j) + 2.0 * f(i - 1, j - 1) -
                  8.0 * f(i, j - 1) + 2.0 * f(i + 1, j - 1) + 1.0 * f(i, j - 2)) / dx4;
  });
}

std::unique_ptr<ManagedArray2D<real_wp>>
Discretization2DCart::createRealArray() const
{
  return std::make_unique<ManagedArray2DOwning<real_wp>>(*this, "temp", ni_, ni_, nhalo_);
}