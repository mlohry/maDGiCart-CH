#include "petsc_time_integrator.hpp"

#include "logger/profiler.hpp"
#include "logger/logger.hpp"

#include <petscvec.h>
#include <iomanip>

//#define PROFILE_PETSC
static std::unique_ptr<Logger::Timer> solve_timer;

namespace {


void
copyPetscVecToSolutionState(Vec x, const IndexArray& indices, SolutionState& state)
{
  profile();

  const PetscScalar* petsc_data;
  VecGetArrayRead(x, &petsc_data);

  auto idx = read_access(indices);
  auto s   = write_access(state.getVec(0));

  maDGForAll(i, 0, idx.size(), {  //
    s[idx[i]] = petsc_data[i];
  });

  VecRestoreArrayRead(x, &petsc_data);
}

void
copySolutionStateToPetscVec(const SolutionState& state, const IndexArray& indices, Vec x)
{
  profile();

  PetscScalar* petsc_data;
  VecGetArray(x, &petsc_data);

  auto idx = read_access(indices);
  auto s   = read_access(state.getVec(0));

  maDGForAll(i, 0, idx.size(), {  //
    petsc_data[i] = s[idx[i]];
  });


  VecRestoreArray(x, &petsc_data);
}


PetscErrorCode
PetscRHSFunction(TS ts, PetscReal time, Vec x, Vec rhs, void* appctx)
{
#ifdef PROFILE_PETSC
  auto timer = Logger::get().timer("PetscRHSFunction");
#endif
  profile();
  PetscTimeIntegrator::Ctx& ctx = *(PetscTimeIntegrator::Ctx*)(appctx);

  copyPetscVecToSolutionState(x, ctx.instance.rhs().interiorIndices(), ctx.instance.getMutableSolutionState());

  ctx.instance.rhs().evalRHSImpl(ctx.instance.getMutableSolutionState(), time, ctx.instance.getMutableResidualState());

  copySolutionStateToPetscVec(ctx.instance.getMutableResidualState(), ctx.instance.rhs().interiorIndices(), rhs);

  return 0;
}


PetscErrorCode
PetscRHSJacobian(TS ts, PetscReal time, Vec x, Mat J, Mat Jpre, void* appctx)
{
#ifdef PROFILE_PETSC
  auto timer = Logger::get().timer("PetscRHSJacobian");
#endif
  profile();

  PetscTimeIntegrator::Ctx& ctx = *(PetscTimeIntegrator::Ctx*)(appctx);

  copyPetscVecToSolutionState(x, ctx.instance.rhs().interiorIndices(), ctx.instance.getMutableSolutionState());

  ctx.instance.rhs().evalJacobian(ctx.instance.getMutableSolutionState(), time, ctx.instance.getJacobian());

  MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  //  Logger::get().TraceMessage("PetscRHSJacobian assembled");
  return 0;
}

PetscErrorCode
PetscTSPreStep(TS ts)
{
#ifdef PROFILE_PETSC
    solve_timer = std::make_unique<Logger::Timer>("Time step", LogLevel::trace);
#endif
  return 0;
}

PetscErrorCode
PetscTSPostStep(TS ts)
{
  const real_t converged_abs_tol = Options::get().converged_abs_tol();

  if (converged_abs_tol) {
    const real_t abs_res = []() {
      const auto& mon = Logger::get().getResidualMonitor();
      for (const auto& kv : mon) {
        if (kv.first == "abs_res") {
          return kv.second;
        }
      }
      Logger::get().FatalMessage("abs_res not found");
      return double(0);
    }();

    std::ostringstream os;
    os << converged_abs_tol;

    if (abs_res <= converged_abs_tol) {
      Logger::get().InfoMessage("Absolute residual convergence " + os.str() + " met.");
      TSSetConvergedReason(ts, TS_CONVERGED_USER);
    }
  }

  const real_t converged_rel_tol = Options::get().converged_rel_tol();

  if (converged_rel_tol) {
    const real_t rel_res = []() {
      const auto& mon = Logger::get().getResidualMonitor();
      for (const auto& kv : mon) {
        if (kv.first == "rel_res") {
          return kv.second;
        }
      }
      Logger::get().FatalMessage("abs_res not found");
      return double(0);
    }();

    std::ostringstream os;
    os << converged_rel_tol;

    if (rel_res <= converged_rel_tol) {
      Logger::get().InfoMessage("Relative residual convergence " + os.str() + " met.");
      TSSetConvergedReason(ts, TS_CONVERGED_USER);
    }
  }


#ifdef PROFILE_PETSC
  solve_timer.reset();
#endif
  return 0;
}


PetscErrorCode
PetscTSMonitor(TS ts, PetscInt step, PetscReal ptime, Vec v, void* appctx)
{
  //  Logger::get().TraceMessage("PetscTSMonitor");
  PetscTimeIntegrator::Ctx& ctx = *(PetscTimeIntegrator::Ctx*)(appctx);
  ctx.instance.notifyObservers(Event::TimeStepComplete);
  ctx.instance.notifyObservers(Event::SolutionUpdate);

  return 0;
}

PetscErrorCode
PetscImplicitSNESMonitor(SNES snes, PetscInt its, PetscReal fnorm, void* appctx)
{
  //  PetscTimeIntegrator::Ctx& ctx = *(PetscTimeIntegrator::Ctx*)(appctx);

  //  Logger::get().TraceMessage("PetscImplicitSNESMonitor");

  return 0;
}

PetscErrorCode
PetscImplicitKSPMonitor(KSP ksp, PetscInt its, PetscReal fnorm, void* appctx)
{
  //  Logger::get().TraceMessage("PetscImplicitKSPMonitor");

  return 0;
}


PetscErrorCode
PetscImplicitLineSearchMonitor(SNESLineSearch ls, void* appctx)
{
  //  Logger::get().TraceMessage("PetscImplicitLineSearchMonitor");

  //  PetscTimeIntegrator::Ctx& ctx = *(PetscTimeIntegrator::Ctx*)(appctx);

  //  Logger::get().TraceMessage(
  //      "Diagonal dominance with time contribution: " +
  //      std::to_string(ctx.instance.getJacobian().diagonalDominance()));

  //  PetscReal norm;
  //  MatNorm(ctx.Jacobian, NORM_1, &norm);
  //  Logger::get().TraceMessage("1 norm: " + std::to_string(norm));

  return 0;
}

}  // namespace


PetscTimeIntegrator::PetscTimeIntegrator(
    TimeIntegrableRHS&           rhs,
    const TimeIntegratorOptions& opts,
    const std::string&           petsc_ts_type)
    : TimeIntegrator(rhs, opts),
      ctx(*this),
      petsc_environment_(PetscEnvironment(Options::get().petsc_options())),
      ts_type_(petsc_ts_type)
{
}

PetscTimeIntegrator::~PetscTimeIntegrator()
{
  if (petscsoln) {
    VecDestroy(&petscsoln);
  }
  if (petscresidual) {
    VecDestroy(&petscresidual);
  }
  if (ctx.ts) {
    TSDestroy(&ctx.ts);
  }
}


void
PetscTimeIntegrator::solve(TimeIntegrableRHS& rhs, InitialConditions& initial_conditions)
{
  profile();

  auto& state0 = getMutableSolutionState();

  initial_conditions.set(rhs, state0);

  initializePetsc(rhs.dofsPerEquation());

  copySolutionStateToPetscVec(state0, ctx.instance.rhs().interiorIndices(), petscsoln);

  TSSetSolution(ctx.ts, petscsoln);
  TSSolve(ctx.ts, petscsoln);

  copyPetscVecToSolutionState(petscsoln, ctx.instance.rhs().interiorIndices(), state0);
}


void
PetscTimeIntegrator::initializePetscJacobian()
{
  TSGetSNES(ctx.ts, &ctx.snes);
  SNESGetKSP(ctx.snes, &ctx.ksp);
  KSPGetPC(ctx.ksp, &ctx.pc);

  auto sparsity_coloring_timer = Logger::get().timer("Creating Jacobian sparsity and matrix coloring");

#ifdef MADG_USE_GPU
  auto      csr   = rhs().createSparseMatrix();
  const int nrows = rhs().dofsPerEquation();
  const int ncols = nrows;

#if defined(MADG_USE_CUDA)
  MatCreateSeqAIJCUSPARSE(PETSC_COMM_WORLD, nrows, ncols, 0, rhs().nNonZerosPerRow().data(), &ctx.Jacobian);
#elif defined(MADG_USE_HIP)
  MatCreateSeqAIJHIPSPARSE(PETSC_COMM_WORLD, nrows, ncols, 0, rhs().nNonZerosPerRow().data(), &ctx.Jacobian);
#endif
  MatSeqAIJSetPreallocationCSR(
      ctx.Jacobian, csr->rowIdx().readHost().data(), csr->colIdx().readHost().data(), csr->values().readHost().data());

  auto         device_borrow_timer = Logger::get().timer("Borrowing CSR pointers on device");
  const int*   d_rowptr;
  const int*   d_colptr;
  real_wp*     d_valptr;
  PetscMemType mtype;
  MatSeqAIJGetCSRAndMemType(ctx.Jacobian, &d_rowptr, &d_colptr, &d_valptr, &mtype);
  printf(
      "Device CSR array pointers from initializePetscJacobian:  %p  %p  %p  %d\n",
      (void*)d_rowptr,
      (void*)d_colptr,
      (void*)d_valptr,
      mtype);

  MatCreate(PETSC_COMM_WORLD, &ctx.JacHost);
  MatSetType(ctx.JacHost, MATSEQAIJ);
  MatSetSizes(ctx.JacHost, nrows, ncols, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSeqAIJSetPreallocationCSR(
      ctx.JacHost, csr->rowIdx().readHost().data(), csr->colIdx().readHost().data(), csr->values().readHost().data());

  const int* h_rowptr;
  const int* h_colptr;
  real_wp*   h_valptr;

  MatSeqAIJGetCSRAndMemType(ctx.JacHost, &h_rowptr, &h_colptr, &h_valptr, &mtype);
  printf(
      "Host CSR array pointers from initializePetscJacobian:  %p  %p  %p\n",
      (void*)h_rowptr,
      (void*)h_colptr,
      (void*)h_valptr);

  csr_jacobian_ = std::make_unique<CSRMatrix>(
      *csr,
      h_valptr,
      const_cast<int*>(h_colptr),
      const_cast<int*>(h_rowptr),
      d_valptr,
      const_cast<int*>(d_colptr),
      const_cast<int*>(d_rowptr));

  device_borrow_timer.elapsed();

#else
  auto      csr   = rhs().createSparseMatrix();
  const int nrows = rhs().dofsPerEquation();
  const int ncols = nrows;

  //  MatCreateSeqAIJ(PETSC_COMM_WORLD, nrows, ncols, 0, rhs().nNonZerosPerRow().data(), &ctx.Jacobian);
  //  MatSeqAIJSetPreallocationCSR(ctx.Jacobian, csr->rowIdx(), csr->colIdx(), csr->values());

  MatCreate(PETSC_COMM_WORLD, &ctx.Jacobian);
  MatSetType(ctx.Jacobian, MATSEQAIJ);
  MatSetSizes(ctx.Jacobian, nrows, ncols, PETSC_DETERMINE, PETSC_DETERMINE);
  MatSeqAIJSetPreallocationCSR(
      ctx.Jacobian, csr->rowIdx().readHost().data(), csr->colIdx().readHost().data(), csr->values().readHost().data());


  auto         host_borrow_timer = Logger::get().timer("Borrowing CSR pointers on host");
  const int*   rowptr;
  const int*   colptr;
  real_wp*     valptr;
  PetscMemType mtype;
  MatSeqAIJGetCSRAndMemType(ctx.Jacobian, &rowptr, &colptr, &valptr, &mtype);
  printf("CSR array pointers from initializePetscJacobian:  %p  %p  %p\n", (void*)rowptr, (void*)colptr, (void*)valptr);

  csr_jacobian_ = std::make_unique<CSRMatrix>(
      *csr, valptr, const_cast<int*>(colptr), const_cast<int*>(rowptr), nullptr, nullptr, nullptr);
  host_borrow_timer.elapsed();
#endif

  //  csr_jacobian_   = rhs().createSparseMatrix();
  //
  //
  //  // initializes a petsc matrix using CSR data.
  //  // this does not copy data -- the petsc matrix is a view into the existing storage
  //
  //  auto matcreatewitharrays = Logger::get().timer("MatCreateSeqAIJWithArrays");
  //  MatCreateSeqAIJWithArrays(
  //      PETSC_COMM_WORLD,
  //      nrows,
  //      ncols,
  //      csr_jacobian_->rowIdx().readHost().data(),
  //      csr_jacobian_->colIdx().readHost().data(),
  //      csr_jacobian_->values().readHost().data(),
  //      &ctx.Jacobian);
  //  matcreatewitharrays.elapsed();


  //  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INFO_DETAIL);
  //  MatView(ctx.Jacobian, PETSC_VIEWER_STDOUT_SELF);


  TSSetRHSJacobian(ctx.ts, ctx.Jacobian, ctx.Jacobian, PetscRHSJacobian, (void*)(&ctx));


  /*
   *
   * run like
   *


   OMP_NUM_THREADS=4 ./maDGiCart --dimension=3 --final_time=10 --max_time_steps=10000 --eps2=0.0019974621629115655
   --sigma=3.5306509073075123 --time_step_size=1e-6 --domain_x_begin=-3.14159 --domain_x_end=3.14159 --log_frequency=1
   --time_integrator=petsc_jacobian_bdf --use_adaptive_time_step=true --time_rel_err_tol=1.e-5
   --petsc_options="-snes_rtol 1e-3 -mg_levels_ksp_max_it 1 -pc_gamg_reuse_interpolation -ts_type bdf -ksp_max_iter 100
   -ksp_atol 1e-6 -ksp_monitor_singular_value  -pc_gamg_agg_nsmooths 0 -mg_levels_ksp_type richardson"

   *
   */

  // for ILU-GMRES
  //  PetscOptionsSetValue(NULL, "-snes_type", "newtonls");
  //  PetscOptionsSetValue(NULL, "-ksp_type", "fgmres");
  //  PetscOptionsSetValue(NULL, "-pc_type", "asm");
  //  PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "100");
  //  PetscOptionsSetValue(NULL, "-sub_pc_factor_levels", "0");


  // for AMG
  PetscOptionsSetValue(NULL, "-snes_type", "newtonls");
  PetscOptionsSetValue(NULL, "-ksp_type", "gmres");
  PetscOptionsSetValue(NULL, "-pc_type", "gamg");
  //  PetscOptionsSetValue(NULL, "-pc_gamg_agg_nsmooths", "0");
  //  PetscOptionsSetValue(NULL, "-mg_levels_ksp_type", "richardson");  // default is chebyshev smoother


  //  PetscOptionsSetValue(NULL, "-ksp_view", NULL);
  //  PetscOptionsSetValue(NULL, "-snes_view", NULL);


  PetscOptionsSetValue(NULL, "-ts_monitor", NULL);
  PetscOptionsSetValue(NULL, "-ksp_monitor", NULL);
  //  PetscOptionsSetValue(NULL, "-ksp_view", NULL);

  PetscOptionsSetValue(NULL, "-ts_adapt_monitor", NULL);
  PetscOptionsSetValue(NULL, "-snes_monitor", NULL);
  PetscOptionsSetValue(NULL, "-snes_converged_reason", NULL);
  PetscOptionsSetValue(NULL, "-ksp_converged_reason", NULL);
  PetscOptionsSetValue(NULL, "-ts_converged_reason", NULL);


  SNESSetFromOptions(ctx.snes);
  KSPSetFromOptions(ctx.ksp);
  PCSetFromOptions(ctx.pc);


  SNESMonitorSet(ctx.snes, PetscImplicitSNESMonitor, &ctx, 0);
  KSPMonitorSet(ctx.ksp, PetscImplicitKSPMonitor, &ctx, 0);
  SNESLineSearch linesearch;
  SNESGetLineSearch(ctx.snes, &linesearch);
  SNESLineSearchMonitorSet(linesearch, PetscImplicitLineSearchMonitor, &ctx, 0);
}


void
PetscTimeIntegrator::initializePetscJFNK()
{
  TSGetSNES(ctx.ts, &ctx.snes);
  SNESGetKSP(ctx.snes, &ctx.ksp);
  KSPGetPC(ctx.ksp, &ctx.pc);

  auto sparsity_coloring_timer = Logger::get().timer("Creating Jacobian sparsity and matrix coloring");

  csr_jacobian_   = rhs().createSparseMatrix();
  const int nrows = rhs().dofsPerEquation();
  const int ncols = nrows;

  // initializes a petsc matrix using CSR data.
  // this does not copy data -- the petsc matrix is a view into the existing storage

  auto matcreatewitharrays = Logger::get().timer("MatCreateSeqAIJWithArrays");
  MatCreateSeqAIJWithArrays(
      PETSC_COMM_WORLD,
      nrows,
      ncols,
      csr_jacobian_->rowIdx().readHost().data(),
      csr_jacobian_->colIdx().readHost().data(),
      csr_jacobian_->values().readHost().data(),
      &ctx.Jacobian);
  matcreatewitharrays.elapsed();

  PetscViewerPushFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INFO_DETAIL);
  MatView(ctx.Jacobian, PETSC_VIEWER_STDOUT_SELF);


  auto        coloring_create = Logger::get().timer("MatColoringCreate");
  MatColoring matcoloring;
  MatColoringCreate(ctx.Jacobian, &matcoloring);
  coloring_create.elapsed();


  MatColoringSetType(matcoloring, MATCOLORINGSL);
  ISColoring iscoloring;
  auto       coloring_apply = Logger::get().timer("MatColoringApply");
  MatColoringApply(matcoloring, &iscoloring);
  coloring_apply.elapsed();

  MatColoringDestroy(&matcoloring);

  sparsity_coloring_timer.elapsed();

  auto          matfdcoloringcreatetimer = Logger::get().timer("Petsc MatFDColoringCreate");
  MatFDColoring fdcoloring;
  MatFDColoringCreate(ctx.Jacobian, iscoloring, &fdcoloring);
  matfdcoloringcreatetimer.elapsed();

  MatFDColoringSetFunction(fdcoloring, (PetscErrorCode(*)(void))SNESTSFormFunction, ctx.ts);

  auto matfdsetuptimer = Logger::get().timer("Petsc MatFDColoringSetUp");
  MatFDColoringSetFromOptions(fdcoloring);
  MatFDColoringSetUp(ctx.Jacobian, iscoloring, fdcoloring);
  matfdsetuptimer.elapsed();


  auto snesetuptimer = Logger::get().timer("Petsc SNESSetJacobian");
  SNESSetJacobian(ctx.snes, ctx.Jacobian, ctx.Jacobian, SNESComputeJacobianDefaultColor, fdcoloring);
  snesetuptimer.elapsed();

  ISColoringDestroy(&iscoloring);
  MatFDColoringDestroy(&fdcoloring);


  PetscOptionsSetValue(NULL, "-snes_mf_operator", NULL);
  PetscOptionsSetValue(NULL, "-snes_type", "newtonls");
  PetscOptionsSetValue(NULL, "-ksp_type", "fgmres");
  PetscOptionsSetValue(NULL, "-pc_type", "asm");
  PetscOptionsSetValue(NULL, "-ksp_gmres_restart", "100");
  PetscOptionsSetValue(NULL, "-sub_pc_factor_levels", "0");

  PetscOptionsSetValue(NULL, "-snes_lag_jacobian", "10000");
  PetscOptionsSetValue(NULL, "-snes_lag_jacobian_persists", "true");
  PetscOptionsSetValue(NULL, "-snes_force_iteration", "true");

  PetscOptionsSetValue(NULL, "-ts_monitor", NULL);
  PetscOptionsSetValue(NULL, "-ksp_monitor", NULL);
  PetscOptionsSetValue(NULL, "-ts_adapt_monitor", NULL);
  PetscOptionsSetValue(NULL, "-snes_monitor", NULL);
  PetscOptionsSetValue(NULL, "-snes_converged_reason", NULL);
  PetscOptionsSetValue(NULL, "-ksp_converged_reason", NULL);
  PetscOptionsSetValue(NULL, "-ts_converged_reason", NULL);


  SNESSetFromOptions(ctx.snes);
  KSPSetFromOptions(ctx.ksp);
  PCSetFromOptions(ctx.pc);
}


void
PetscTimeIntegrator::initializePetsc(int ndofs)
{
  auto initialize_petsc_timer = Logger::get().timer("PetscTimeIntegrator initialization");

#ifdef MADG_USE_GPU
#ifdef MADG_USE_CUDA
  auto cuda_malloc_timer = Logger::get().timer("Allocating Petsc CUDA arrays");
  VecCreateSeqCUDA(PETSC_COMM_WORLD, ndofs, &petscsoln);
  VecCreateSeqCUDA(PETSC_COMM_WORLD, ndofs, &petscresidual);
  cuda_malloc_timer.elapsed();
#endif
#ifdef MADG_USE_HIP
  auto hip_malloc_timer = Logger::get().timer("Allocating Petsc HIP arrays");
  VecCreateSeqHIP(PETSC_COMM_WORLD, ndofs, &petscsoln);
  VecCreateSeqHIP(PETSC_COMM_WORLD, ndofs, &petscresidual);
  hip_malloc_timer.elapsed();
#endif
#else
  VecCreateSeq(PETSC_COMM_WORLD, ndofs, &petscsoln);
  VecCreateSeq(PETSC_COMM_WORLD, ndofs, &petscresidual);
#endif

  TSCreate(PETSC_COMM_WORLD, &ctx.ts);
  TSSetProblemType(ctx.ts, TS_NONLINEAR);
  TSSetRHSFunction(ctx.ts, NULL, PetscRHSFunction, &ctx);

  if (ts_type_ == "ode23") {
    Logger::get().InfoMessage("Initializing petsc with ode23 integrator");
    TSSetType(ctx.ts, TSRK);
    TSRKSetType(ctx.ts, TSRK3BS);
  }
  else if (ts_type_ == "jfnk_bdf") {
    Logger::get().InfoMessage("Initializing petsc with JFNK BDF integrator");
    initializePetscJFNK();
    TSSetType(ctx.ts, TSBDF);
  }
  else if (ts_type_ == "jacobian_bdf") {
    Logger::get().InfoMessage("Initializing petsc with explicit Jacobian BDF integrator");
    initializePetscJacobian();
    //    TSSetType(ctx.ts, TSBDF);
  }
  else {
    Logger::get().WarningMessage("Petsc TS type " + ts_type_ + " not recognized.");
  }

  const auto& time_opts = getTimeOptions();
  TSSetTimeStep(ctx.ts, PetscReal(time_opts.dt_initial_));
  TSSetTime(ctx.ts, PetscReal(time_opts.t0_));
  TSSetStepNumber(ctx.ts, time_opts.initial_step_);


  if (time_opts.use_adaptive_time_step_) {
    TSAdapt adapt;
    TSGetAdapt(ctx.ts, &adapt);
    TSAdaptSetType(adapt, TSADAPTBASIC);
    TSAdaptSetStepLimits(adapt, time_opts.min_time_step_size_, time_opts.max_time_step_size_);
    TSSetTolerances(ctx.ts, time_opts.time_rel_err_tol_, NULL, time_opts.time_rel_err_tol_, NULL);
  }
  else {
    TSAdapt adapt;
    TSGetAdapt(ctx.ts, &adapt);
    TSAdaptSetType(adapt, TSADAPTNONE);
    TSSetTimeStep(ctx.ts, PetscReal(Options::get().time_step_size()));
  }

  if (time_opts.max_time_steps_ != 0) {
    TSSetMaxSteps(ctx.ts, time_opts.max_time_steps_);
  }
  else {
    TSSetMaxSteps(ctx.ts, std::numeric_limits<PetscInt>::max());
  }

  if (time_opts.tfinal_ != 0) {
    TSSetMaxTime(ctx.ts, time_opts.tfinal_);
  }
  else {
    TSSetMaxTime(ctx.ts, std::numeric_limits<PetscReal>::max());
  }

  TSSetExactFinalTime(ctx.ts, TS_EXACTFINALTIME_MATCHSTEP);

  TSSetPostStep(ctx.ts, PetscTSPostStep);
  TSSetPreStep(ctx.ts, PetscTSPreStep);
  TSMonitorSet(ctx.ts, PetscTSMonitor, &ctx, 0);

  TSSetFromOptions(ctx.ts);
}


real_t
PetscTimeIntegrator::getCurrentTime() const
{
  PetscReal t;
  TSGetTime(ctx.ts, &t);
  return (double)t;
}

real_t
PetscTimeIntegrator::getTimeStepSize() const
{
  PetscReal dt;
  TSGetTimeStep(ctx.ts, &dt);
  return (double)dt;
}

int_t
PetscTimeIntegrator::getCurrentStep() const
{
  PetscInt iter;
  TSGetStepNumber(ctx.ts, &iter);
  return (int_t)iter;
}


static auto petsc_ode23instance = FactoryRegistry<TimeIntegrator>::get().add(
    "petsc_ode23",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<PetscTimeIntegrator>(rhs, opts, "ode23");
    });

static auto petsc_jfnk_bdf_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "petsc_jfnk_bdf",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<PetscTimeIntegrator>(rhs, opts, "jfnk_bdf");
    });

static auto petsc_jacobian_bdf_instance = FactoryRegistry<TimeIntegrator>::get().add(
    "petsc_jacobian_bdf",
    [](TimeIntegrableRHS& rhs, const TimeIntegratorOptions& opts) {
      return std::make_unique<PetscTimeIntegrator>(rhs, opts, "jacobian_bdf");
    });
