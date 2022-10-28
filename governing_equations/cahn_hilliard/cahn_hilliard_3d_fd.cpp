#include "cahn_hilliard_3d_fd.hpp"

#include "logger/logger.hpp"
#include "program_options/program_options.hpp"


CahnHilliard3DFD::CahnHilliard3DFD(Discretization3DCart& geom, const CahnHilliardParameters& params)
    : CahnHilliardBase(params), geom_(geom), use_fused_kernels_(Options::get().use_fused_kernels())
{
  laplacian_rhs_term_ = geom_.createRealArray();
  laplacian_argument_ = geom_.createRealArray();
  biharmonic_term_    = geom_.createRealArray();
  linear_term_        = geom_.createRealArray();
}


void
CahnHilliard3DFD::evalRHSImpl(const SolutionState& flovars, double time, SolutionState& rhs)
{
  if (use_fused_kernels_){
    evalRHSFused(flovars, time, rhs);
  } else {
    evalRHSBasic(flovars, time, rhs);
  }
}


void
CahnHilliard3DFD::evalRHSBasic(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const CahnHilliardState3D& state = dynamic_cast<const CahnHilliardState3D&>(flovars);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  geom_.applyNeumannBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  /**
   * Compute  laplacian(u*c^3 -b*c)
   */
  {


    auto c     = read_access(state.c().asArray());
    auto del2f = write_access(laplacian_argument_->asArray());

    maDGForAll(i, 0, c.size(), { del2f[i] = pow(c[i], 3.0) - c[i]; });

    geom_.laplacian(*laplacian_argument_, *laplacian_rhs_term_);
  }

  /**
   * Compute biharmonic(c)
   */
  geom_.biharmonic(state.c(), *biharmonic_term_);

  /**
   * compute sigma*(c-m)
   */
  {
    auto          c     = read_access(state.c().asArray());
    auto          idx   = read_access(geom_.interiorIndices());
    auto          term  = write_access(linear_term_->asArray());
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

    auto rhseval = write_access(dstate_dt.c().asArray());
    auto idx     = read_access(geom_.interiorIndices());
    auto del4    = read_access(biharmonic_term_->asArray());
    auto del2    = read_access(laplacian_rhs_term_->asArray());
    auto sigterm = read_access(linear_term_->asArray());

    const real_wp eps2 = this->eps2();

    maDGForAll(ii, 0, idx.size(), {
      const int i = idx[ii];
      rhseval[i]  = -eps2 * del4[i] + del2[i] - sigterm[i];
    });
  }
}

void
CahnHilliard3DFD::evalRHSFused(const SolutionState& flovars, double time, SolutionState& rhs)
{
  profile();
  const CahnHilliardState3D& state     = dynamic_cast<const CahnHilliardState3D&>(flovars);
  CahnHilliardState3D&       dstate_dt = dynamic_cast<CahnHilliardState3D&>(rhs);

  geom_.applyPeriodicBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));
  geom_.applyNeumannBoundaryConditions(const_cast<ManagedArray3DNonOwning<real_wp>&>(state.c()));

  /**
   * Compute  laplacian(u*c^3 -b*c)
   */
  {
    auto c     = read_access(state.c().asArray());
    auto del2f = write_access(laplacian_argument_->asArray());
    maDGForAll(i, 0, c.size(), { del2f[i] = pow(c[i], 3.0) - c[i]; });
  }

  {
    auto          laplacian_arg = read_access((*laplacian_argument_));
    auto          idx           = read_access(interiorIndices());
    const real_wp dx2           = geom_.dx() * geom_.dx();
    const real_wp dx4           = dx2 * dx2;
    auto          f             = read_access(state.c());

    auto rhseval = write_access(dstate_dt.c());

    const real_wp m     = this->m();
    const real_wp sigma = this->sigma();
    const real_wp eps2  = this->eps2();

    maDGForAll(ii, 0, idx.size(), {
      int i;
      int j;
      int k;
      laplacian_arg.getIJK(idx[ii], i, j, k);


      const real_wp del2 = (laplacian_arg(i - 1, j, k) + laplacian_arg(i + 1, j, k) + laplacian_arg(i, j - 1, k) +
                            laplacian_arg(i, j + 1, k) + laplacian_arg(i, j, k - 1) + laplacian_arg(i, j, k + 1) -
                            real_wp(6) * laplacian_arg(i, j, k)) /
                           dx2;

      const real_wp del4f =
          (real_wp(42) * f(i, j, k) -
           real_wp(12) *
               (f(i + 1, j, k) + f(i, j + 1, k) + f(i, j, k + 1) + f(i - 1, j, k) + f(i, j - 1, k) + f(i, j, k - 1)) +
           f(i + 2, j, k) + f(i, j + 2, k) + f(i, j, k + 2) + f(i - 2, j, k) + f(i, j - 2, k) + f(i, j, k - 2) +
           real_wp(2) * (f(i - 1, j - 1, k) + f(i - 1, j + 1, k) + f(i + 1, j - 1, k) + f(i + 1, j + 1, k) +
                         f(i - 1, j, k - 1) + f(i + 1, j, k - 1) + f(i, j - 1, k - 1) + f(i, j + 1, k - 1) +
                         f(i - 1, j, k + 1) + f(i + 1, j, k + 1) + f(i, j - 1, k + 1) + f(i, j + 1, k + 1))) /
          dx4;

      const real_wp linear_term = sigma * (f(i, j, k) - m);

      rhseval(i, j, k) = -eps2 * del4f + del2 - linear_term;
    });
  }
}
