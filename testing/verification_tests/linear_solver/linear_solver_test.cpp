#include "data_structures/csr_matrix.hpp"

#include <gtest/gtest.h>

#include "data_structures/ldu_csr_matrix.hpp"
#include "file_io/read_csr_matrix.hpp"

TEST(LinearSolver, LDUCSRJacobi)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto spmat  = readCSRMatrix(row_file, col_file, val_file);
  auto ldumat = std::make_unique<LDUCSRMatrix>(*spmat);

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->jacobiSmooth(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}

TEST(LinearSolver, CSRJacobi)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto spmat = readCSRMatrix(row_file, col_file, val_file);

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = spmat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    spmat->jacobiSmooth(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = spmat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}


TEST(LinearSolver, CSRJacobiCacheDiagonal)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto spmat = readCSRMatrix(row_file, col_file, val_file);

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = spmat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    spmat->jacobiSmoothCacheDiagonalInverse(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = spmat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}

TEST(LinearSolver, CSRJacobiBLASCacheDiagonal)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto spmat = readCSRMatrix(row_file, col_file, val_file);

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = spmat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    spmat->jacobiSmoothBLASCacheDiagonalInverse(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = spmat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}

TEST(LinearSolver, LDUCSRCSRJacobiCacheDiagonal)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

//  auto spmat = readCSRMatrix(row_file, col_file, val_file);
  auto ldumat = std::make_unique<LDUCSRMatrix>(*readCSRMatrix(row_file, col_file, val_file));

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->jacobiSmoothCacheDiagonalInverse(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}

TEST(LinearSolver, LDUCSRJacobiBLASCacheDiagonal)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  //  auto spmat = readCSRMatrix(row_file, col_file, val_file);
  auto ldumat = std::make_unique<LDUCSRMatrix>(*readCSRMatrix(row_file, col_file, val_file));

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 0.25;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->jacobiSmoothBLASCacheDiagonalInverse(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNorm(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 103200, 500);
}

TEST(LinearSolver, LDUCSRForwardGaussSeidelBLAS)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto ldumat = std::make_unique<LDUCSRMatrix>(*readCSRMatrix(row_file, col_file, val_file));

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 1.0;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->forwardGaussSeidelSmooth(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNormDiagScaledRHS(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 13501, 500);
}


TEST(LinearSolver, LDUCSRBackwardGaussSeidelBLAS)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto ldumat = std::make_unique<LDUCSRMatrix>(*readCSRMatrix(row_file, col_file, val_file));

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 1.0;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->backwardGaussSeidelSmooth(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNormDiagScaledRHS(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 13501, 500);
}


TEST(LinearSolver, LDUCSRSymmetricGaussSeidelBLAS)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

  const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto ldumat = std::make_unique<LDUCSRMatrix>(*readCSRMatrix(row_file, col_file, val_file));

  ManagedArray<real_wp> sol_array({}, "test_sol", rhs_array.size(), real_wp(0));

  const double abs_err_0 = ldumat->rmsNorm(sol_array, rhs_array);
  double       rel_err   = 1.0;
  std::cout << "err: " << abs_err_0 << "  rel_err: " << rel_err << "\n";

  const double weight      = 1.0;
  int          iter        = 0;
  const int    max_iter    = 200000;
  const double rel_err_tol = 1.e-14;

  do {
    ldumat->symmetricGaussSeidelSmooth(sol_array, rhs_array, weight);
    if (iter % 100 == 0) {
      const double abs_err = ldumat->rmsNormDiagScaledRHS(sol_array, rhs_array);
      rel_err              = abs_err / abs_err_0;
      std::cout << "iter: " << iter << " err: " << abs_err << "  rel_err: " << rel_err << "\n";
    }
    iter++;
  } while (rel_err > rel_err_tol && iter < max_iter);

  EXPECT_NEAR(rel_err, 1.e-14, 1.e-14);
  EXPECT_NEAR(iter, 10501, 500);
}
