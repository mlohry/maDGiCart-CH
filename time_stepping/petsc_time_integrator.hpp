#pragma once

#include <petscksp.h>
#include <petscsnes.h>
#include <petscts.h>

#include "petsc_environment.hpp"
#include "time_integrator.hpp"


class PetscTimeIntegrator : public TimeIntegrator {
 public:
  PetscTimeIntegrator(TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts, const std::string& petsc_ts_type);

  ~PetscTimeIntegrator();

  void solve(TimeIntegrableRHS&, InitialConditions&) override;
  void doTimeStep(TimeIntegrableRHS&, SolutionState&, SolutionState&, double time, double dt) override
  { /* no-op */
  }

  int_t  getCurrentStep() const override;
  real_t getCurrentTime() const override;
  real_t getTimeStepSize() const override;

  SolutionState& getMutableSolutionState()
  {
    return const_cast<SolutionState&>(TimeIntegrator::getCurrentSolutionState());
  }
  SolutionState& getMutableResidualState() { return const_cast<SolutionState&>(TimeIntegrator::getCurrentResidual()); }

  struct Ctx {
    Ctx(PetscTimeIntegrator& outer) : instance(outer) {}

    TS   ts       = nullptr;
    Mat  Jacobian = nullptr;
    Mat  JacHost  = nullptr;
    SNES snes     = nullptr;
    KSP  ksp      = nullptr;
    PC   pc       = nullptr;

    PetscTimeIntegrator& instance;
  };

  CSRMatrix& getJacobian() const { return *csr_jacobian_; }

 private:
  Ctx                    ctx;
  const PetscEnvironment petsc_environment_;
  const std::string      ts_type_;
  const std::string      ts_specific_type_;

  Vec petscsoln     = nullptr;
  Vec petscresidual = nullptr;

  std::unique_ptr<CSRMatrix> csr_jacobian_;


  void initializePetsc(int ndofs);
  void initializePetscJFNK();
  void initializePetscJacobian();
};
