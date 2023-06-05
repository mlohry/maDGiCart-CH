#include "celero/Celero.h"
#include "data_structures/csr_matrix.hpp"
#include "data_structures/ldu_csr_matrix.hpp"
#include "file_io/read_csr_matrix.hpp"


class SparseMatrixFixture : public celero::TestFixture {
 public:
  SparseMatrixFixture()
  {
    const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
    const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
    const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
    const std::string rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";

    const std::vector<real_wp> rhs_vec = readArray<real_wp>(rhs_file);
    rhs_array_                         = std::make_unique<ManagedArray<real_wp>>(
        ManagedArrayOwner{}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

    csr_mat_     = readCSRMatrix(row_file, col_file, val_file);
    ldu_csr_mat_ = std::make_unique<LDUCSRMatrix>(*csr_mat_);
    sol_array_ =
        std::make_unique<ManagedArray<real_wp>>(ManagedArrayOwner{}, "test_sol", rhs_array_->size(), real_wp(0));
  }


 protected:
  std::unique_ptr<LDUCSRMatrix>          ldu_csr_mat_;
  std::unique_ptr<CSRMatrix>             csr_mat_;
  std::unique_ptr<ManagedArray<real_wp>> rhs_array_;
  std::unique_ptr<ManagedArray<real_wp>> sol_array_;
  const real_wp                          weight_ = 0.25;
};

static constexpr int SamplesCount    = 1;  // set to 0 for automatic determination
static constexpr int IterationsCount = 500;


BASELINE_F(JacobiSmooth, CSR, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  csr_mat_->jacobiSmooth(*sol_array_, *rhs_array_, weight_);
//  csr_mat_->jacobiSmoothCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}

BENCHMARK_F(JacobiSmooth, CSRCacheDiag, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  csr_mat_->jacobiSmoothCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}

BENCHMARK_F(JacobiSmooth, CSRBlas, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  csr_mat_->jacobiSmoothBLASCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}

BENCHMARK_F(JacobiSmooth, LDUCSR, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  ldu_csr_mat_->jacobiSmooth(*sol_array_, *rhs_array_, weight_);
//  ldu_csr_mat_->jacobiSmoothCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}

BENCHMARK_F(JacobiSmooth, LDUCSRCacheDiag, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  ldu_csr_mat_->jacobiSmoothCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}

BENCHMARK_F(JacobiSmooth, LDUCSRBlas, SparseMatrixFixture, SamplesCount, IterationsCount)
{  //
  ldu_csr_mat_->jacobiSmoothBLASCacheDiagonalInverse(*sol_array_, *rhs_array_, weight_);
}
