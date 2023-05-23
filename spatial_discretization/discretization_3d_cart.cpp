#include "discretization_3d_cart.hpp"

#include "data_structures/scalar_solution_state.hpp"
#include "data_structures/solution_state.hpp"
#include "program_options/program_options.hpp"


Discretization3DCart::Discretization3DCart(const CartesianDomainDefinition& domain)
    : SpatialDiscretization("Discretization3DCart", domain),
      ni_(domain.nx),
      nj_(domain.ny),
      nk_(domain.nz),
      nhalo_(domain.nhalo),
      ninhalo_(ni_ + 2 * nhalo_),
      njnhalo_(nj_ + 2 * nhalo_),
      nknhalo_(nk_ + 2 * nhalo_),
      dx_((domain.xend - domain.xbeg) / real_wp(ni_)),
      interior_indices_(*this, "InteriorIndices", ni_ * nj_ * nk_),
      periodic_donor_indices_(*this, "PeriodicDonorIndices", ninhalo_ * njnhalo_ * nknhalo_ - ni_ * nj_ * nk_),
      periodic_receiver_indices_(*this, "PeriodicReceiverIndices", ninhalo_ * njnhalo_ * nknhalo_ - ni_ * nj_ * nk_),
      xmin_indices_(*this, "XMinIndices", njnhalo_ * nknhalo_ * nhalo_),
      xmax_indices_(*this, "XMaxIndices", njnhalo_ * nknhalo_ * nhalo_),
      ymin_indices_(*this, "YMinIndices", ninhalo_ * nknhalo_ * nhalo_),
      ymax_indices_(*this, "YMaxIndices", ninhalo_ * nknhalo_ * nhalo_),
      zmin_indices_(*this, "ZMinIndices", ninhalo_ * njnhalo_ * nhalo_),
      zmax_indices_(*this, "ZMaxIndices", ninhalo_ * njnhalo_ * nhalo_),
      x_coord_(*this, "CoordinateX", ni_, nj_, nk_, nhalo_),
      y_coord_(*this, "CoordinateY", ni_, nj_, nk_, nhalo_),
      z_coord_(*this, "CoordinateZ", ni_, nj_, nk_, nhalo_),
      x_vertex_coord_(*this, "VertexCoordinateX", ni_ + 1, nj_ + 1, nk_ + 1, nhalo_),
      y_vertex_coord_(*this, "VertexCoordinateY", ni_ + 1, nj_ + 1, nk_ + 1, nhalo_),
      z_vertex_coord_(*this, "VertexCoordinateZ", ni_ + 1, nj_ + 1, nk_ + 1, nhalo_)
{
  {
    auto idx_list = write_access_host(interior_indices_);
    int  inode    = 0;

    for (int i = 0; i < ni_; ++i) {
      for (int j = 0; j < nj_; ++j) {
        for (int k = 0; k < nk_; ++k) {
          idx_list[inode] = get1DindexFrom3D(i, j, k, nhalo_, njnhalo_, nknhalo_);
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
      for (int j = -nhalo_; j < nj_ + nhalo_ + 1; ++j) {
        for (int k = -nhalo_; k < nk_ + nhalo_ + 1; ++k) {
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
      for (int j = -nhalo_; j < nj_ + nhalo_; ++j) {
        for (int k = -nhalo_; k < nk_ + nhalo_; ++k) {
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
      for (int j = -nhalo_; j < nj_ + nhalo_; ++j) {
        for (int k = -nhalo_; k < nk_ + nhalo_; ++k) {
          bool is_halo_cell = [&]() {
            if (i < 0 || i >= ni_ || j < 0 || j >= nj_ || k < 0 || k >= nk_) {
              return true;
            }
            return false;
          }();

          if (is_halo_cell) {
            receiver_list[inode] = get1DindexFrom3D(i, j, k, nhalo_, njnhalo_, nknhalo_);

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
                return j + nj_;
              }
              if (j >= nj_) {
                return j - nj_;
              }
              return j;
            }();

            const int kdonor = [&]() {
              if (k < 0) {
                return k + nk_;
              }
              if (k >= nk_) {
                return k - nk_;
              }
              return k;
            }();

            donor_list[inode] = get1DindexFrom3D(idonor, jdonor, kdonor, nhalo_, njnhalo_, nknhalo_);

            inode++;
          }
        }
      }
    }
  }
}


std::unique_ptr<SpatialDiscretization>
Discretization3DCart::createCoarsenedDiscretization() const
{
  const auto fine_domain = this->domain();
  auto       domain      = fine_domain;
  if ((domain.nx % 2) || (domain.ny % 2) || (domain.nz % 2)) {
    Logger::get().FatalMessage("createCoarsenedDiscretization needs domain size divisible by 2");
  }

  domain.nx /= 2;
  domain.ny /= 2;
  domain.nz /= 2;
  return std::make_unique<Discretization3DCart>(domain);
}


void
Discretization3DCart::interpolateFineToCoarse(
    const SolutionState&         fine,
    const SpatialDiscretization& coarse_geom,
    SolutionState&               coarse) const
{
  profile();
  auto coarse_idx   = read_access(dynamic_cast<const Discretization3DCart&>(coarse_geom).interiorIndices());
  auto fine_state   = read_access(dynamic_cast<const ScalarSolutionState3D&>(fine).c());
  auto coarse_state = write_access(dynamic_cast<ScalarSolutionState3D&>(coarse).c());

  maDGForAll(ii, 0, coarse_idx.size(), {
    int i;
    int j;
    int k;
    coarse_state.getIJK(coarse_idx[ii], i, j, k);

    const int         finei  = 2 * i;
    const int         finej  = 2 * j;
    const int         finek  = 2 * k;
    constexpr real_wp factor = 1. / 8.;

    coarse_state(i, j, k) =
        factor * (fine_state(finei, finej, finek) + fine_state(finei + 1, finej, finek) +
                  fine_state(finei + 1, finej + 1, finek) + fine_state(finei, finej + 1, finek) +
                  fine_state(finei, finej, finek + 1) + fine_state(finei + 1, finej, finek + 1) +
                  fine_state(finei + 1, finej + 1, finek + 1) + fine_state(finei, finej + 1, finek + 1));
  });
}


void
Discretization3DCart::interpolateCoarseToFine(
    const SolutionState&         coarse,
    const SpatialDiscretization& fine_geom,
    SolutionState&               fine) const
{
  profile();
  auto coarse_idx   = read_access(this->interiorIndices());
  auto coarse_state = read_access(dynamic_cast<const ScalarSolutionState3D&>(coarse).c());
  auto fine_state   = write_access(dynamic_cast<ScalarSolutionState3D&>(fine).c());

  if (Options::get().mg_interpolation() == "nearest") {  // 0th order interpolation, probably shouldn't be used

    maDGForAll(ii, 0, coarse_idx.size(), {
      int i;
      int j;
      int k;
      coarse_state.getIJK(coarse_idx[ii], i, j, k);

      const int finei = 2 * i;
      const int finej = 2 * j;
      const int finek = 2 * k;

      const real_wp val = coarse_state(i, j, k);

      fine_state(finei, finej, finek)             = val;
      fine_state(finei + 1, finej, finek)         = val;
      fine_state(finei + 1, finej + 1, finek)     = val;
      fine_state(finei, finej + 1, finek)         = val;
      fine_state(finei, finej, finek + 1)         = val;
      fine_state(finei + 1, finej, finek + 1)     = val;
      fine_state(finei + 1, finej + 1, finek + 1) = val;
      fine_state(finei, finej + 1, finek + 1)     = val;
    });
  }
  else if (Options::get().mg_interpolation() == "linear") {  // trilinear interpolation

    const int ni = this->ni();
    const int nj = this->nj();
    const int nk = this->nk();

    maDGForAll(i, -1, ni, {
      const int finei = 2 * i + 1;
      for (int j = -1; j < nj; ++j) {
        const int finej = 2 * j + 1;
        for (int k = -1; k < nk; ++k) {
          const int finek = 2 * k + 1;

          const real_wp s_ijk   = coarse_state(i, j, k);
          const real_wp s_i1jk  = coarse_state(i + 1, j, k);
          const real_wp s_i1j1k = coarse_state(i + 1, j + 1, k);
          const real_wp s_ij1k  = coarse_state(i, j + 1, k);

          const real_wp s_ijk1   = coarse_state(i, j, k + 1);
          const real_wp s_i1jk1  = coarse_state(i + 1, j, k + 1);
          const real_wp s_i1j1k1 = coarse_state(i + 1, j + 1, k + 1);
          const real_wp s_ij1k1  = coarse_state(i, j + 1, k + 1);

          fine_state(finei, finej, finek) = 0.140625 * s_i1j1k + 0.046875 * s_i1j1k1 + 0.140625 * s_i1jk +
                                            0.046875 * s_i1jk1 + 0.046875 * s_ij1k + 0.015625 * s_ij1k1 +
                                            0.421875 * s_ijk + 0.140625 * s_ijk1;

          fine_state(finei + 1, finej, finek) = 0.046875 * s_i1j1k + 0.015625 * s_i1j1k1 + 0.421875 * s_i1jk +
                                                0.140625 * s_i1jk1 + 0.140625 * s_ij1k + 0.046875 * s_ij1k1 +
                                                0.140625 * s_ijk + 0.046875 * s_ijk1;

          fine_state(finei + 1, finej + 1, finek) = 0.140625 * s_i1j1k + 0.046875 * s_i1j1k1 + 0.140625 * s_i1jk +
                                                    0.046875 * s_i1jk1 + 0.421875 * s_ij1k + 0.140625 * s_ij1k1 +
                                                    0.046875 * s_ijk + 0.015625 * s_ijk1;

          fine_state(finei, finej + 1, finek) = 0.421875 * s_i1j1k + 0.140625 * s_i1j1k1 + 0.046875 * s_i1jk +
                                                0.015625 * s_i1jk1 + 0.140625 * s_ij1k + 0.046875 * s_ij1k1 +
                                                0.140625 * s_ijk + 0.046875 * s_ijk1;

          fine_state(finei, finej, finek + 1) = 0.046875 * s_i1j1k + 0.140625 * s_i1j1k1 + 0.046875 * s_i1jk +
                                                0.140625 * s_i1jk1 + 0.015625 * s_ij1k + 0.046875 * s_ij1k1 +
                                                0.140625 * s_ijk + 0.421875 * s_ijk1;

          fine_state(finei + 1, finej, finek + 1) = 0.015625 * s_i1j1k + 0.046875 * s_i1j1k1 + 0.140625 * s_i1jk +
                                                    0.421875 * s_i1jk1 + 0.046875 * s_ij1k + 0.140625 * s_ij1k1 +
                                                    0.046875 * s_ijk + 0.140625 * s_ijk1;

          fine_state(finei + 1, finej + 1, finek + 1) = 0.046875 * s_i1j1k + 0.140625 * s_i1j1k1 + 0.046875 * s_i1jk +
                                                        0.140625 * s_i1jk1 + 0.140625 * s_ij1k + 0.421875 * s_ij1k1 +
                                                        0.015625 * s_ijk + 0.046875 * s_ijk1;

          fine_state(finei, finej + 1, finek + 1) = 0.140625 * s_i1j1k + 0.421875 * s_i1j1k1 + 0.015625 * s_i1jk +
                                                    0.046875 * s_i1jk1 + 0.046875 * s_ij1k + 0.140625 * s_ij1k1 +
                                                    0.046875 * s_ijk + 0.140625 * s_ijk1;
        }
      }
    });
  }
  else {
    Logger::get().FatalMessage("mg_interpolation not recognized");
  }
}


std::unique_ptr<ManagedArray3D<real_wp>>
Discretization3DCart::createRealArray() const
{
  return std::make_unique<ManagedArray3DOwning<real_wp>>(*this, "temp", ni_, nj_, nk_, nhalo_);
}


void
Discretization3DCart::laplacian(const ManagedArray3D<real_wp>& state_in, ManagedArray3D<real_wp>& del2state_out) const
{
  profile();
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
  if (domain().xbc == BCType::Neumann) {
    applyNeumannBCX(state);
  }
  if (domain().ybc == BCType::Neumann) {
    applyNeumannBCY(state);
  }
  if (domain().zbc == BCType::Neumann) {
    applyNeumannBCZ(state);
  }
}

void
Discretization3DCart::applyNeumannBCX(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto nhalo   = nhalo_;
  const auto ni      = ni_;
  const auto nj      = nj_;
  const auto nk      = nk_;

  maDGForAll(j, -nhalo, nj + nhalo, {
    for (int k = -nhalo; k < nk + nhalo; ++k) {
      for (int i = -nhalo; i < 0; ++i) {
        state(i, j, k) = state(0, j, k);
      }
      for (int i = ni; i < ni + nhalo; ++i) {
        state(i, j, k) = state(ni - 1, j, k);
      }
    }
  });
}

void
Discretization3DCart::applyNeumannBCY(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto nhalo   = nhalo_;
  const auto ni      = ni_;
  const auto nj      = nj_;
  const auto nk      = nk_;

  maDGForAll(i, -nhalo, ni + nhalo, {
    for (int k = -nhalo; k < nk + nhalo; ++k) {
      for (int j = -nhalo; j < 0; ++j) {
        state(i, j, k) = state(i, 0, k);
      }
      for (int j = nj; j < nj + nhalo; ++j) {
        state(i, j, k) = state(i, nj - 1, k);
      }
    }
  });
}

void
Discretization3DCart::applyNeumannBCZ(ManagedArray3D<real_wp>& state_in) const
{
  profile();
  auto state = read_write_access(state_in);

  const auto nhalo   = nhalo_;
  const auto ni      = ni_;
  const auto nj      = nj_;
  const auto nk      = nk_;

  maDGForAll(i, -nhalo, ni + nhalo, {
    for (int j = -nhalo; j < nj + nhalo; ++j) {
      for (int k = -nhalo; k < 0; ++k) {
        state(i, j, k) = state(i, j, 0);
      }
      for (int k = nk; k < nk + nhalo; ++k) {
        state(i, j, k) = state(i, j, nk - 1);
      }
    }
  });
}
