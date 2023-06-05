#pragma once

#include "ode23_smoother.hpp"
#include "spatial_discretization/spatial_discretization.hpp"
#include "time_integrator.hpp"


class BDFMultigrid : public TimeIntegrator {
 public:
  BDFMultigrid(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts);

  void doTimeStep(TimeIntegrableRHS& rhs, SolutionState& state, SolutionState& dstate_dt, double time, double dt)
      override;

  int    nCyclesToConverge() const { return n_cycles_to_converge_; }
  double convergenceRate() const { return convergence_rate_; }

  class FineGridRHS : public TimeIntegrableRHS {
   public:
    FineGridRHS(TimeIntegrableRHS& rhs, BDFMultigrid& parent);

    void setSolOld(const SolutionState& sol_old);

    void evalRHSImpl(const SolutionState& state, double time, SolutionState& dstate_dt) override;

    void evalJacobian(const SolutionState& state, double time, CSRMatrix& J) override
    {
      Logger::get().FatalMessage("not implemented");
    }

    std::vector<std::string>           equationNames() const override { return rhs_.equationNames(); }
    int_t                              nEquations() const override { return rhs_.nEquations(); }
    int_t                              dofsPerEquation() const override { return rhs_.dofsPerEquation(); }
    std::unique_ptr<SolutionState>     createSolutionState() const override { return rhs_.createSolutionState(); }
    std::unique_ptr<CSRMatrix>         createSparseMatrix() const override { return rhs_.createSparseMatrix(); }
    std::vector<int>                   nNonZerosPerRow() const override { return rhs_.nNonZerosPerRow(); }
    const IndexArray&                  interiorIndices() const override { return rhs_.interiorIndices(); }
    const SpatialDiscretization*       hasGeometry() const override { return rhs_.hasGeometry(); }
    double                             cflScale() const override { return rhs_.cflScale(); }
    std::unique_ptr<TimeIntegrableRHS> clone(SpatialDiscretization& g) const override { return rhs_.clone(g); }

   private:
    TimeIntegrableRHS&             rhs_;
    SolutionState const*           sol_old_;
    std::unique_ptr<SolutionState> res0_;
    BDFMultigrid&                  parent_;
  };

  class CoarsenedLevelRHS : public TimeIntegrableRHS {
   public:
    CoarsenedLevelRHS(const TimeIntegrableRHS& fine_rhs, const SpatialDiscretization& fine_geom, BDFMultigrid& parent);

    void evalRHSImpl(const SolutionState& state, double time, SolutionState& dstate_dt) override;

    void evalJacobian(const SolutionState& state, double time, CSRMatrix& J) override
    {
      Logger::get().FatalMessage("not implemented");
    }

    std::vector<std::string>       equationNames() const override { return rhs_->equationNames(); }
    int_t                          nEquations() const override { return rhs_->nEquations(); }
    int_t                          dofsPerEquation() const override { return rhs_->dofsPerEquation(); }
    std::unique_ptr<SolutionState> createSolutionState() const override { return rhs_->createSolutionState(); }
    std::unique_ptr<CSRMatrix>     createSparseMatrix() const override { return rhs_->createSparseMatrix(); }
    std::vector<int>               nNonZerosPerRow() const override { return rhs_->nNonZerosPerRow(); }
    const IndexArray&              interiorIndices() const override { return rhs_->interiorIndices(); }
    const SpatialDiscretization*   hasGeometry() const override { return rhs_->hasGeometry(); }
    double                         cflScale() const override { return rhs_->cflScale(); }

    void setCoarseStatesFromFine(const SolutionState& fine_state, const SolutionState& fine_residual, double time);

    void applyRelaxedDelta(real_wp relax, const SolutionState& delta, SolutionState& state);
    void getUpdateDeltaOnFine(SolutionState& fine_delta);

    std::unique_ptr<SolutionState> state_;
    std::unique_ptr<SolutionState> state_presmooth_;
    std::unique_ptr<SolutionState> residual_presmooth_;
    std::unique_ptr<SolutionState> restricted_fine_grid_residual_;
    std::unique_ptr<SolutionState> dstate_dt_;
    std::unique_ptr<SolutionState> update_delta_;
    std::unique_ptr<SolutionState> mg_source_;
    std::unique_ptr<SolutionState> fine_grid_delta_;  // try to remove this

    std::unique_ptr<TimeIntegrableRHS> clone(SpatialDiscretization& g) const override { return fine_rhs_.clone(g); }

    //   private:
    const TimeIntegrableRHS&     fine_rhs_;  // only used for clone
    const SpatialDiscretization& fine_geom_;
    BDFMultigrid&                parent_;

    std::unique_ptr<TimeIntegrableRHS>     rhs_;
    std::unique_ptr<SpatialDiscretization> geom_;
  };

  real_wp physicalTimeStep() const { return dt_; }


 private:
  std::unique_ptr<FineGridRHS>   backward_euler_rhs_;
  std::unique_ptr<SolutionState> sol_old_;

  std::vector<TimeIntegrableRHS*> mg_rhs_;
  std::vector<std::unique_ptr<ODE23Smoother>> mg_smoothers_;
  std::vector<std::unique_ptr<CoarsenedLevelRHS>> coarse_rhs_vec_;

  real_wp dt_;
  real_wp time_;

  const double           rel_tol_;
  const int              max_cycles_;
  const int              verbosity_;
  const std::vector<int> mg_cycle_definition_;

  int    n_cycles_to_converge_ = 0;
  double convergence_rate_     = 0;

  void doTimeStepCycle(
      const std::vector<int>& cycle,
      TimeIntegrableRHS&      rhs,
      SolutionState&          state,
      SolutionState&          dstate_dt,
      double                  time,
      double                  dt);

  static std::vector<int> getMGCycle();
};
