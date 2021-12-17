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
CahnHilliard2DFD::setRandomInitialCondition(SolutionState& state) const
{
  std::random_device rd;  // Will be used to obtain a seed for the random number engine
  //  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::mt19937                            gen(1);  // Standard mersenne_twister_engine seeded with 1
  std::uniform_real_distribution<real_wp> dist(initial_value_min_, initial_value_max_);


  auto idx = read_access_host(geom_.interiorIndices());
  auto c   = write_access_host(state.getVec(0));

  /**
   * Would prefer to use a maDGForAllHost here but std::uniform_real_distribution lambda copy
   * has to be mutable.
   */
  for (int i = 0; i < idx.size(); ++i) {
    c[idx[i]] = dist(gen);
  }
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


std::map<std::string, double>
CahnHilliard2DFD::solutionReport(const SolutionState& state, const SolutionState& residual)
{
  std::map<std::string, double> reports;

  {
    ReduceSumReal squared_norm(0);
    auto          rhs = read_access(residual.getVec(0));
    auto          idx = read_access(geom_.interiorIndices());

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      squared_norm += pow(rhs[i], 2.0);
    });
    reports["residual"] = std::sqrt(squared_norm.get() / geom_.nInteriorPoints());
  }

  //  {
  //    ReduceSumReal gradient_magnitude(0);
  //    auto f = read_access(dynamic_cast<const CahnHilliardState&>(state).c());
  //    auto idx = read_access(geom_.interiorIndices());
  //    const real_wp dx = geom_.dx();
  //    const real_wp dy = geom_.dy();
  //
  //    maDGForAll(ii, 0, idx.size(), {
  //      int i;
  //      int j;
  //      f.getIJ(idx[ii], i, j);
  //      wrong expression, this is a laplacian, dummy
  //      const real_wp grad = (f(i-1,j) - 2.0*f(i,j) + f(i+1,j))/(dx*dx) + (f(i,j+1) - 2.0*f(i,j) + f(i,j-1))/(dy*dy);
  //      gradient_magnitude += abs(grad)*dx*dy;
  //    });
  //
  //    reports["gradient_magnitude"] = gradient_magnitude.get();
  //  }


  return reports;
}

