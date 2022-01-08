#include "cahn_hilliard_2d_fd.hpp"

#include <random>


CahnHilliard2DFD::CahnHilliard2DFD(
    Discretization2DCart&         geom,
    const CahnHilliardParameters& params)
    : geom_(geom),
      m_(params.m()),
      sigma_(params.sigma()),
      eps2_(params.eps2()),
      initial_value_min_(params.initialMin()),
      initial_value_max_(params.initialMax())
{
}


void
CahnHilliard2DFD::evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const CahnHilliardState& state = dynamic_cast<const CahnHilliardState&>(flovars);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray2DNonOwning<real_wp>&>(state.c()));

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
    const real_wp m     = m_;
    const real_wp sigma = sigma_;

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      term[i]     = sigma * (c[i] - m);
    });
  }


  {
    CahnHilliardState& dstate_dt = dynamic_cast<CahnHilliardState&>(rhs);
    {
      auto r = write_access(dstate_dt.c().asArray());
      maDGForAll(i, 0, r.size(), { r[i] = real_wp(0); });
    }

    auto rhs     = write_access(dstate_dt.c().asArray());
    auto idx     = read_access(geom_.interiorIndices());
    auto del4    = read_access(biharmonic_term->asArray());
    auto del2    = read_access(laplacian_rhs_term->asArray());
    auto sigterm = read_access(linear_term->asArray());

    const real_wp eps2 = eps2_;


    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      rhs[i]      = -eps2 * del4[i] + del2[i] - sigterm[i];
    });
  }
}
