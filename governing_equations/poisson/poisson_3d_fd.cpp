#include "poisson_3d_fd.hpp"


Poisson3DFD::Poisson3DFD(Discretization3DCart& geom) : geom_(geom) {}


void
Poisson3DFD::evalRHSImpl(const SolutionState& state_in, double time, SolutionState& rhs)
{
  profile();
  ScalarSolutionState3D& state =
      const_cast<ScalarSolutionState3D&>(dynamic_cast<const ScalarSolutionState3D&>(state_in));
  ScalarSolutionState3D& dstate_dt = dynamic_cast<ScalarSolutionState3D&>(rhs);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  this->applyDirichletBoundaryConditions(state);

  geom_.laplacian(state.c(), dstate_dt.c());
  addSourceTerm(dstate_dt);
}

void
Poisson3DFD::addSourceTerm(ScalarSolutionState3D& rhs_in)
{
  profile();
  auto x_arr = read_access(geom_.x());
  auto y_arr = read_access(geom_.y());
  auto z_arr = read_access(geom_.z());
  auto rhs   = read_write_access(rhs_in.c());
  auto idx   = read_access(geom_.interiorIndices());

  maDGForAll(ii, 0, idx.size(), {
    int i;
    int j;
    int k;
    rhs.getIJK(idx[ii], i, j, k);

    const real_t x = x_arr(i, j, k);
    const real_t y = y_arr(i, j, k);
    const real_t z = z_arr(i, j, k);

    rhs(i, j, k) += 3.0 * sin(x) * cos(y) * cos(z);
  });
}

void
Poisson3DFD::applyDirichletBoundaryConditions(ScalarSolutionState3D& state_in)
{
  profile();

  auto state = read_write_access(state_in.c());
  auto x_arr = read_access(geom_.x());
  auto y_arr = read_access(geom_.y());
  auto z_arr = read_access(geom_.z());

  const auto nhalo   = geom_.nhalo();
  const auto nk      = geom_.nk();
  const auto ni      = geom_.ni();
  const auto nj      = geom_.nj();

  maDGForAll(i, -nhalo, ni + nhalo, {
    for (int j = -nhalo; j < nj + nhalo; ++j) {
      for (int k = -nhalo; k < 0; ++k) {
        const real_t x = x_arr(i, j, k);
        const real_t y = y_arr(i, j, k);
        const real_t z = z_arr(i, j, k);
        state(i, j, k) = sin(x) * cos(y) * cos(z);
      }
      for (int k = nk; k < nk + nhalo; ++k) {
        const real_t x = x_arr(i, j, k);
        const real_t y = y_arr(i, j, k);
        const real_t z = z_arr(i, j, k);
        state(i, j, k) = sin(x) * cos(y) * cos(z);
      }
    }
  })
}


std::unique_ptr<TimeIntegrableRHS>
Poisson3DFD::clone(SpatialDiscretization& geom) const
{
  return std::make_unique<Poisson3DFD>(dynamic_cast<Discretization3DCart&>(geom));
}
