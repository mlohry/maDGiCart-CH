#include "cahn_hilliard_3d_fd.hpp"

#include "logger/logger.hpp"

CahnHilliard3DFD::CahnHilliard3DFD(Discretization3DCart& geom, const CahnHilliardParameters& params)
: CahnHilliardBase(params), geom_(geom)
{
}


void CahnHilliard3DFD::evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const CahnHilliardState3D& state = dynamic_cast<const CahnHilliardState3D&>(flovars);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  /**
   * Compute  laplacian(u*c^3 -b*c)
   */
  auto laplacian_rhs_term = geom_.createRealArray();
  {
    auto laplacian_argument = geom_.createRealArray();

    auto c     = read_access(state.c().asArray());
    auto del2f = write_access(laplacian_argument->asArray());

    maDGForAll(i, 0, c.size(), { del2f[i] = pow(c[i], 3.0) - c[i]; });

    geom_.laplacian(*laplacian_argument, *laplacian_rhs_term);
  }

  /**
   * Compute biharmonic(c)
   */
  auto biharmonic_term = geom_.createRealArray();
  geom_.biharmonic(state.c(), *biharmonic_term);

  /**
   * compute sigma*(c-m)
   */
  auto linear_term = geom_.createRealArray();
  {
    auto          c     = read_access(state.c().asArray());
    auto          idx   = read_access(geom_.interiorIndices());
    auto          term  = write_access(linear_term->asArray());
    const real_wp m     = this->m();
    const real_wp sigma = this->sigma();

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      term[i]     = sigma * (c[i] - m);
    });
  }

  {
    CahnHilliardState3D& dstate_dt = dynamic_cast<CahnHilliardState3D&>(rhs);
    {
      auto r = write_access(dstate_dt.c().asArray());
      maDGForAll(i, 0, r.size(), { r[i] = real_wp(0); });
    }

    auto rhs     = write_access(dstate_dt.c().asArray());
    auto idx     = read_access(geom_.interiorIndices());
    auto del4    = read_access(biharmonic_term->asArray());
    auto del2    = read_access(laplacian_rhs_term->asArray());
    auto sigterm = read_access(linear_term->asArray());

    const real_wp eps2 = this->eps2();


    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      rhs[i]      = -eps2 * del4[i] + del2[i] - sigterm[i];
    });
  }
}
