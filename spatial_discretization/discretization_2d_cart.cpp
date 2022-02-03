#include "discretization_2d_cart.hpp"

#include <memory>

Discretization2DCart::Discretization2DCart(const CartesianDomainDefinition& domain)
    : SpatialDiscretization("Discretization2DCart", domain),
      ni_(domain.nx),
      nhalo_(domain.nhalo),
      ninhalo_(ni_ + 2 * nhalo_),
      dx_((domain.xend - domain.xbeg) / real_wp(ni_)),
      interior_indices_(*this, "InteriorIndices", ni_ * ni_),
      boundary_indices_(*this, "BoundaryIndices", (ni_ + 2 * nhalo_) * nhalo_ * 2 + ni_ * nhalo_ * 2),
      x_coord_(*this, "CoordinateX", ni_, ni_, nhalo_),
      y_coord_(*this, "CoordinateY", ni_, ni_, nhalo_),
      x_vertex_coord_(*this, "VertexCoordinateX", ni_ + 1, ni_ + 1, nhalo_),
      y_vertex_coord_(*this, "VertexCoordinateY", ni_ + 1, ni_ + 1, nhalo_)
{
  {
    auto idx_list = write_access_host(interior_indices_);
    int  inode    = 0;

    for (int i = 0; i < ni_; ++i) {
      for (int j = 0; j < ni_; ++j) {
        idx_list[inode] = get1Dindex(i, j, nhalo_, ninhalo_);
        inode++;
      }
    }
  }

  {
    int ibcnode = 0;
    // bottom rows including corners
    auto bc_list = write_access_host(boundary_indices_);
    for (int j = -nhalo_; j < 0; ++j) {
      for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // top rows including corners
    for (int j = ni_; j < ni_ + nhalo_; ++j) {
      for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // left excluding corners
    for (int i = -nhalo_; i < 0; ++i) {
      for (int j = 0; j < ni_; ++j) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }

    // right excluding corners
    for (int i = ni_; i < ni_ + nhalo_; ++i) {
      for (int j = 0; j < ni_; ++j) {
        bc_list[ibcnode] = get1Dindex(i, j, nhalo_, ninhalo_);
        ibcnode++;
      }
    }
  }

  {
    auto xv = write_access_host(x_vertex_coord_);
    auto yv = write_access_host(y_vertex_coord_);

    for (int i = -nhalo_; i < ni_ + nhalo_ + 1; ++i) {
      for (int j = -nhalo_; j < ni_ + nhalo_ + 1; ++j) {
        xv(i, j) = domain.xbeg + dx_ * real_wp(i);
        yv(i, j) = domain.ybeg + dx_ * real_wp(j);
      }
    }
  }

  {
    auto xv = read_access_host(x_vertex_coord_);
    auto yv = read_access_host(y_vertex_coord_);
    auto x  = write_access_host(x_coord_);
    auto y  = write_access_host(y_coord_);

    for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
      for (int j = -nhalo_; j < ni_ + nhalo_; ++j) {
        x(i, j) = 0.5 * (xv(i, j) + xv(i + 1, j));
        y(i, j) = 0.5 * (yv(i, j) + yv(i, j + 1));
      }
    }
  }
}


void
Discretization2DCart::printState(const ManagedArray2D<real_wp>& state_in) const
{
  auto state = read_access_host(state_in);

  for (int j = ni_ + nhalo_ - 1; j >= -nhalo_; --j) {
    for (int i = -nhalo_; i < ni_ + nhalo_; ++i) {
      std::cout << std::scientific << state(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

void
Discretization2DCart::applyPeriodicBoundaryConditions(ManagedArray2D<real_wp>& state_in) const
{
  profile();
  auto      state  = read_write_access(state_in);
  const int nhalos = nhalo();
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
                   8.0 * f(i, j - 1) + 2.0 * f(i + 1, j - 1) + 1.0 * f(i, j - 2)) /
                  dx4;
  });
}

std::unique_ptr<ManagedArray2D<real_wp>>
Discretization2DCart::createRealArray() const
{
  return std::make_unique<ManagedArray2DOwning<real_wp>>(*this, "temp", ni_, ni_, nhalo_);
}