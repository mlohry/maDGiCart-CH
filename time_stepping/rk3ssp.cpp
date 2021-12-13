#include "rk3ssp.hpp"

#include "data_structures/exec_includes.hpp"


RK3SSP::RK3SSP(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) : TimeIntegrator(rhs, opts)
{
  k1_rhs_         = rhs.createSolutionState();
  k2_rhs_         = rhs.createSolutionState();
  k3_rhs_         = rhs.createSolutionState();
  stage_solution_ = rhs.createSolutionState();
}


void RK3SSP::doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
{
  profile();

  rhs.evalRHSImpl(state, time, *k1_rhs_);

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto stage_sol = write_access(stage_solution_->getVec(ivec).asArray());
    auto sol0      = read_access(state.getVec(ivec).asArray());
    auto k1        = read_access(k1_rhs_->getVec(ivec).asArray());

    const real_wp dt  = dt;
    const real_wp a21 = a21_;

    maDGForAll(i, 0, stage_sol.size(), {  //
      stage_sol[i] = sol0[i] + dt * a21 * k1[i];
    });
  }

  rhs.evalRHSImpl(*stage_solution_, time + c2_ * dt, *k2_rhs_);

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto stage_sol = write_access(stage_solution_->getVec(ivec).asArray());
    auto sol0      = read_access(state.getVec(ivec).asArray());
    auto k1        = read_access(k1_rhs_->getVec(ivec).asArray());
    auto k2        = read_access(k2_rhs_->getVec(ivec).asArray());

    const real_wp dt  = dt;
    const real_wp a31 = a31_;
    const real_wp a32 = a32_;

    maDGForAll(i, 0, stage_sol.size(), {  //
      stage_sol[i] = sol0[i] + dt * (a31 * k1[i] + a32 * k2[i]);
    });
  }

  rhs.evalRHSImpl(*stage_solution_, time + c3_ * dt, *k3_rhs_);

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto residual = write_access(dstate_dt.getVec(ivec).asArray());
    auto k1       = read_access(k1_rhs_->getVec(ivec).asArray());
    auto k2       = read_access(k2_rhs_->getVec(ivec).asArray());
    auto k3       = read_access(k3_rhs_->getVec(ivec).asArray());

    const real_wp b1 = b1_;
    const real_wp b2 = b2_;
    const real_wp b3 = b3_;

    maDGForAll(i, 0, residual.size(), {  //
      residual[i] = b1 * k1[i] + b2 * k2[i] + b3 * k3[i];
    });
  }

  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
    auto sol = read_write_access(state.getVec(ivec).asArray());
    auto res = read_access(dstate_dt.getVec(ivec).asArray());

    maDGForAll(i, 0, sol.size(), {  //
      sol[i] += dt * res[i];
    });
  }
}


static auto rk3sspinstance = FactoryRegistry<TimeIntegrator>::get().add(
    "rk3ssp",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) { return std::make_unique<RK3SSP>(rhs, opts); });
