#include "bdf_multigrid.hpp"

#include "ode23_smoother.hpp"

#include <iomanip>


BDFMultigrid::BDFMultigrid(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts)
    : TimeIntegrator(rhs, opts),
      rel_tol_(Options::get().mg_rel_tol()),
      max_cycles_(Options::get().mg_max_cycles()),
      verbosity_(Options::get().mg_verbosity()),
      mg_cycle_definition_(getMGCycle())
{
  backward_euler_rhs_ = std::make_unique<FineGridRHS>(rhs, *this);
  mg_smoothers_.push_back(std::make_unique<ODE23Smoother>(*backward_euler_rhs_));
  mg_rhs_.push_back(backward_euler_rhs_.get());

  for (int ilevel = 2; ilevel <= Options::get().mg_levels(); ++ilevel) {
    coarse_rhs_vec_.push_back(
        std::make_unique<CoarsenedLevelRHS>(*mg_rhs_.back(), *mg_rhs_.back()->hasGeometry(), *this));
    mg_rhs_.push_back(coarse_rhs_vec_.back().get());
    mg_smoothers_.push_back(std::make_unique<ODE23Smoother>(*coarse_rhs_vec_.back()));
  }
}


namespace {
void
copySolution(const SolutionState& in, SolutionState& out)
{
  for (int ivec = 0; ivec < in.nvecs(); ++ivec) {
    auto vec_out = write_access(out.getVec(ivec));
    auto vec_in  = read_access(in.getVec(ivec));
    maDGForAll(i, 0, vec_in.size(), { vec_out[i] = vec_in[i]; });
  }
}

std::vector<double>
reduce_residual(const TimeIntegrableRHS& rhs, const SolutionState& residuals)
{
  std::vector<double> res(residuals.nvecs());

  auto idx = read_access(rhs.interiorIndices());

  for (int ieqn = 0; ieqn < residuals.nvecs(); ++ieqn) {
    ReduceSumReal squared_norm(0);
    auto          array = read_access(residuals.getVec(ieqn));
    maDGForAll(ii, 0, idx.size(), { squared_norm += pow(array[idx[ii]], 2.0); });
    res[ieqn] = sqrt(squared_norm.get() / idx.size());
  }

  return res;
}
}  // namespace


void
BDFMultigrid::doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
{
  time_    = time;
  dt_      = dt;
  sol_old_ = rhs.createSolutionState();
  copySolution(state, *sol_old_);
  backward_euler_rhs_->setSolOld(*sol_old_);

  doTimeStepCycle(mg_cycle_definition_, rhs, state, dstate_dt, time, dt);
}

void
BDFMultigrid::doTimeStepCycle(
    const std::vector<int>& cycle,
    TimeIntegrableRHS&      rhs,
    SolutionState&          state,
    SolutionState&          dstate_dt,
    double                  time,
    double                  dt)
{
  // todo we are wasting a residual evaluation here
  backward_euler_rhs_->evalRHSImpl(state, time, dstate_dt);
  const auto res0 = reduce_residual(*backward_euler_rhs_, dstate_dt);
  auto       resN = res0;

  const auto start_clock = std::chrono::high_resolution_clock::now();

  if (verbosity_) {
    std::stringstream ss;
    ss << std::setw(8) << "iter" << std::setw(15) << "residual";
    for (int ilevel = 0; ilevel < mg_smoothers_.size(); ++ilevel) {
      ss << std::setw(12) << "dt lvl " + std::to_string(ilevel);
    }
    ss << "\n";
    ss << "-----------------------------------------------------------------------\n";
    ss << std::setw(8) << 0;
    ss << std::setw(15) << std::scientific << std::setprecision(5) << res0[0];

    for (const auto& mg_smoother : mg_smoothers_) {
      ss << std::setw(12) << std::scientific << std::setprecision(3) << mg_smoother->stepSize();
    }
    Logger::get().TraceMessage(ss.str());
  }

  const real_wp coarse_update_relaxation = Options::get().mg_coarse_relaxation();
  const real_wp fine_update_relaxation   = Options::get().mg_relaxation();

  for (int icycle = 0; icycle < max_cycles_; ++icycle) {

    int ilevel = 0;

    for (const auto& mg_stage : cycle) {

      if (mg_stage == 0) {  // smoothing stage
        if (ilevel == 0) {
          mg_smoothers_[ilevel]->doSmooth(*mg_rhs_[ilevel], state, dstate_dt, time);
        }
        else {
          mg_smoothers_[ilevel]->doSmooth(
              *mg_rhs_[ilevel],
              *dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel])->state_,
              *dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel])->dstate_dt_,
              time);
        }
      }
      else if (mg_stage == 1) {  // transfer fine to coarse
        if (ilevel == 0) {
          dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel + 1])
              ->setCoarseStatesFromFine(state, mg_smoothers_[ilevel]->getLastResidual(), time);
        }
        else {
          dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel + 1])
              ->setCoarseStatesFromFine(
                  *dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel])->state_,
                  mg_smoothers_[ilevel]->getLastResidual(),
                  time);
        }
        ilevel++;
      }
      else if (mg_stage == -1) {  // transfer coarse to fine
        auto* lvl_rhs = dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel]);
        lvl_rhs->getUpdateDeltaOnFine(*lvl_rhs->fine_grid_delta_);

        if (ilevel - 1 == 0) {
          lvl_rhs->applyRelaxedDelta(fine_update_relaxation, *lvl_rhs->fine_grid_delta_, state);
        }
        else {
          lvl_rhs->applyRelaxedDelta(
              coarse_update_relaxation,
              *lvl_rhs->fine_grid_delta_,
              *dynamic_cast<CoarsenedLevelRHS*>(mg_rhs_[ilevel - 1])->state_);
        }
        ilevel--;
      }
    }

    resN = reduce_residual(*backward_euler_rhs_, dstate_dt);

    if (verbosity_) {
      std::stringstream ss;
      ss << std::setw(8) << icycle + 1;
      ss << std::setw(15) << std::scientific << std::setprecision(5) << resN[0];
      for (const auto& mg_smoother : mg_smoothers_) {
        ss << std::setw(12) << std::scientific << std::setprecision(3) << mg_smoother->stepSize();
      }
      Logger::get().TraceMessage(ss.str());
    }

    const double convergence = resN[0] / res0[0];

    if (convergence <= rel_tol_) {
      n_cycles_to_converge_ = icycle + 1;
      convergence_rate_     = std::pow(convergence, 1.0 / double(n_cycles_to_converge_));
      std::stringstream ss;
      if (verbosity_) {
        ss << "-----------------------------------------------------------------------\n";
      }
      std::chrono::duration<double> dt   = std::chrono::high_resolution_clock::now() - start_clock;
      ss << "Multigrid converged to relative tolerance.\n"
         << "Relative error:   " << std::scientific << std::setprecision(4) << convergence << "\n"
         << "Number of cycles: " << n_cycles_to_converge_ << "\n"
         << "Convergence rate: " << std::scientific << std::setprecision(4) << convergence_rate_ << "\n"
         << "Time (s):         " << std::scientific << std::setprecision(4) << dt.count();
      Logger::get().TraceMessage(ss.str());
      return;
    }
  }

  std::stringstream ss;
  if (verbosity_) {
    ss << "-----------------------------------------------------------------------\n";
  }
  ss << "Multigrid failed to converged in " << max_cycles_ << " cycles. Relative residual " << std::scientific
     << resN[0] / res0[0];
  n_cycles_to_converge_ = 0;
  convergence_rate_     = std::pow(resN[0] / res0[0], 1.0 / double(max_cycles_));
  Logger::get().WarningMessage(ss.str());
}


BDFMultigrid::FineGridRHS::FineGridRHS(TimeIntegrableRHS& rhs, BDFMultigrid& parent)
    : rhs_(rhs), res0_(createSolutionState()), parent_(parent)
{
}


void
BDFMultigrid::FineGridRHS::setSolOld(const SolutionState& sol_old)
{
  sol_old_ = &sol_old;
  rhs_.evalRHSImpl(sol_old, parent_.time_, *res0_);
}


void
BDFMultigrid::FineGridRHS::evalRHSImpl(const SolutionState& state, double time, SolutionState& dstate_dt)
{
  rhs_.evalRHSImpl(state, time, dstate_dt);

  if (Options::get().time_integrator() != "steady") {

    for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
      auto          res_n = read_write_access(dstate_dt.getVec(ivec));
      auto          sol0  = read_access(sol_old_->getVec(ivec));
      auto          sol   = read_access(state.getVec(ivec));
      const real_wp dt    = parent_.dt_;

      maDGForAll(i, 0, res_n.size(), {  //
        res_n[i] = res_n[i] - (1.0 / dt) * (sol[i] - sol0[i]);
      });
    }
  }
}


BDFMultigrid::CoarsenedLevelRHS::CoarsenedLevelRHS(
    const TimeIntegrableRHS&     fine_rhs,
    const SpatialDiscretization& fine_geom,
    BDFMultigrid&                parent)
    : fine_rhs_(fine_rhs), fine_geom_(fine_geom), parent_(parent)
{
  geom_ = fine_geom_.createCoarsenedDiscretization();
  rhs_  = fine_rhs.clone(*geom_);

  state_                         = rhs_->createSolutionState();
  state_                         = rhs_->createSolutionState();
  state_presmooth_               = rhs_->createSolutionState();
  dstate_dt_                     = rhs_->createSolutionState();
  update_delta_                  = rhs_->createSolutionState();
  mg_source_                     = rhs_->createSolutionState();
  restricted_fine_grid_residual_ = rhs_->createSolutionState();
  residual_presmooth_            = rhs_->createSolutionState();
  fine_grid_delta_               = fine_rhs.createSolutionState();
}

void
BDFMultigrid::CoarsenedLevelRHS::evalRHSImpl(const SolutionState& state, double time, SolutionState& dstate_dt)
{
  rhs_->evalRHSImpl(state, time, dstate_dt);

  if (mg_source_) {
    for (int ivec = 0; ivec < state.nvecs(); ++ivec) {
      auto parent_residual  = read_write_access(dstate_dt.getVec(ivec));
      auto multigrid_source = read_access(mg_source_->getVec(ivec));
      maDGForAll(i, 0, parent_residual.size(), { parent_residual[i] += multigrid_source[i]; });
    }
  }
}

void
BDFMultigrid::CoarsenedLevelRHS::setCoarseStatesFromFine(
    const SolutionState& fine_state,
    const SolutionState& fine_residual,
    double               time)
{
  fine_geom_.interpolateFineToCoarse(fine_state, *geom_, *state_);
  fine_geom_.interpolateFineToCoarse(fine_residual, *geom_, *restricted_fine_grid_residual_);

  // store the initial coarse grid solution to compute the delta afterwards
  for (int ivec = 0; ivec < state_->nvecs(); ++ivec) {
    auto s  = read_access(state_->getVec(ivec));
    auto s0 = write_access(state_presmooth_->getVec(ivec));
    maDGForAll(i, 0, s.size(), { s0[i] = s[i]; });
  }

  // rhs_evaluator_ evaluates the spatial residual without source terms
  rhs_->evalRHSImpl(*state_, time, *dstate_dt_);  // time missing?

  // compute the coarse grid multigrid source term
  for (int ivec = 0; ivec < dstate_dt_->nvecs(); ++ivec) {
    auto coarse_src = write_access((mg_source_->getVec(ivec)));
    auto fine_res   = read_access((restricted_fine_grid_residual_->getVec(ivec)));
    auto coarse_res = read_access(dstate_dt_->getVec(ivec));

    maDGForAll(i, 0, coarse_src.size(), { coarse_src[i] = -coarse_res[i] + fine_res[i]; })
  }

  this->evalRHSImpl(*state_presmooth_, time, *residual_presmooth_);

  // todo probably needs unsteady correction
}


void
BDFMultigrid::CoarsenedLevelRHS::applyRelaxedDelta(real_wp relax, const SolutionState& delta, SolutionState& state)
{
  for (int ivec = 0; ivec < state.nvecs(); ++ivec) {

    auto update = read_access(delta.getVec(ivec));
    auto soln   = read_write_access(state.getVec(ivec));

    const real_wp r = relax;

    maDGForAll(i, 0, update.size(), { soln[i] += r * update[i]; });
  }
}

void
BDFMultigrid::CoarsenedLevelRHS::getUpdateDeltaOnFine(SolutionState& fine_delta)
{
  for (int ivec = 0; ivec < state_->nvecs(); ++ivec) {

    auto coarse0 = read_access(state_presmooth_->getVec(ivec));
    auto coarse  = read_access(state_->getVec(ivec));
    auto delta   = write_access(update_delta_->getVec(ivec));

    maDGForAll(i, 0, coarse.size(), { delta[i] = coarse[i] - coarse0[i]; });
  }

  geom_->interpolateCoarseToFine(*update_delta_, fine_geom_, fine_delta);
}


std::vector<int>
BDFMultigrid::getMGCycle()
{
  std::vector<int>       cycle;
  const int              nlevels = Options::get().mg_levels();
  const int              nsmooth = Options::get().mg_nsmooth();
  const std::vector<int> smoothing(nsmooth, 0);

  switch (nlevels) {
    case 1:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;
    case 2:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;
    case 3:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;
    case 4:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;
    case 5:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;
    case 6:
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      cycle.push_back(-1);
      cycle.insert(cycle.end(), smoothing.begin(), smoothing.end());
      break;

    default:
      Logger::get().FatalMessage("mg_levels " + std::to_string(Options::get().mg_levels()) + " not supported.");
  }

  return cycle;
}

static auto bdfinstance =
    FactoryRegistry<TimeIntegrator>::get().add("bdf1", [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<BDFMultigrid>(rhs, opts);
    });

static auto steadyinstance =
    FactoryRegistry<TimeIntegrator>::get().add("steady", [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<BDFMultigrid>(rhs, opts);
    });
