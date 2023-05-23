#include "data_structures/csr_matrix.hpp"

#include <gtest/gtest.h>

#include "data_structures/ldu_csr_matrix.hpp"
#include "file_io/read_csr_matrix.hpp"


namespace {
std::unique_ptr<CSRMatrix>
construct1DPoissonMatrix(int nrows)
{
  const int ncols = nrows;
  const int nnz   = nrows * 3 - 4;  // -4 because we will make the first and last rows only have 1 on the diagonal
  auto      spmat = std::make_unique<CSRMatrix>(nrows, ncols, nnz);

  // set up the row indexes, which is of length nrows+1
  {
    auto row_index = write_access_host(spmat->rowIdx());

    maDGForAllHost(irow, 0, nrows + 1, {
      // first row_index is always 0
      if (irow == 0) {
        row_index[irow] = 0;
      }
      else if (irow == nrows) {
        row_index[irow] = nnz;
      }
      else {
        const int nnz_above_this_row = irow * 3 - 2;
        row_index[irow]              = nnz_above_this_row;
      }
    });
  }

  // set up the col indexes and values
  {
    auto col_index = write_access_host(spmat->colIdx());
    auto row_index = read_access_host(spmat->rowIdx());
    auto vals      = write_access_host(spmat->values());


    maDGForAllHost(irow, 0, nrows, {
      const int row_start = row_index[irow];
      const int row_end   = row_index[irow + 1];
      if (irow == 0) {
        col_index[row_start] = 0;
        vals[row_start]      = 1.0;
      }
      else if (irow == nrows - 1) {
        col_index[row_start] = nrows - 1;
        vals[row_start]      = 1.0;
      }
      else {
        ASSERT_EQ(row_end, row_start + 3);
        col_index[row_start]     = irow - 1;
        col_index[row_start + 1] = irow;
        col_index[row_start + 2] = irow + 1;
        vals[row_start]          = -1.0;
        vals[row_start + 1]      = 2.0;
        vals[row_start + 2]      = -1.0;
      }
    });
  }

  return spmat;
}


std::unique_ptr<CSRMatrix>
construct5x5TestMatrix()
{
  /*
   * constuct the matrix
   * [1, 2, 0, 11, 0]
   * [0, 3, 4, 0, 0]
   * [0, 5, 6, 7, 0]
   * [0, 0, 0, 8, 0]
   * [0, 0, 0, 9, 10]
   */
  auto spmat = std::make_unique<CSRMatrix>(5, 5, 11);

  auto col_index = write_access_host(spmat->colIdx());
  auto row_index = write_access_host(spmat->rowIdx());
  auto vals      = write_access_host(spmat->values());

  row_index[0] = 0;
  row_index[1] = 3;
  row_index[2] = 5;
  row_index[3] = 8;
  row_index[4] = 9;
  row_index[5] = 11;

  col_index[0]  = 0;
  col_index[1]  = 1;
  col_index[2]  = 3;
  col_index[3]  = 1;
  col_index[4]  = 2;
  col_index[5]  = 1;
  col_index[6]  = 2;
  col_index[7]  = 3;
  col_index[8]  = 3;
  col_index[9]  = 3;
  col_index[10] = 4;

  vals[0]  = 1;
  vals[1]  = 2;
  vals[2]  = 11;
  vals[3]  = 3;
  vals[4]  = 4;
  vals[5]  = 5;
  vals[6]  = 6;
  vals[7]  = 7;
  vals[8]  = 8;
  vals[9]  = 9;
  vals[10] = 10;

  return spmat;
}

}  // namespace

TEST(CSRMatrix, SparseMatrixVectorMultiplication)
{
  // test correctness for the 1d poisson matrix
  {
    const int nrows = 10;
    auto      spmat = construct1DPoissonMatrix(nrows);

    ManagedArray<real_wp> x({}, "x", nrows, 1.0);
    ManagedArray<real_wp> b({}, "b", nrows);


    spmat->matVecMultiply(x, b);

    auto result = read_access_host(b);
    maDGForAllHost(irow, 0, nrows, {
      if (irow == 0) {
        ASSERT_EQ(result[irow], 1.0) << " at row : " << irow;
      }
      else if (irow == nrows - 1) {
        ASSERT_EQ(result[irow], 1.0) << " at row : " << irow;
      }
      else {
        ASSERT_EQ(result[irow], 0.0) << " at row : " << irow;
      }
    });
  }

  // test correctness for the 5x5 test matrix
  {
    auto      spmat = construct5x5TestMatrix();
    const int nrows = spmat->rows();

    ManagedArray<real_wp> x({}, "x", nrows);
    ManagedArray<real_wp> b({}, "b", nrows);

    auto x_arr = write_access_host(x);
    maDGForAllHost(irow, 0, nrows, { x_arr[irow] = irow + 1; });
    // set x={1,2,3,4,5}

    spmat->matVecMultiply(x, b);

    auto result = read_access_host(b);
    ASSERT_EQ(result[0], 49.0);
    ASSERT_EQ(result[1], 18.0);
    ASSERT_EQ(result[2], 56.0);
    ASSERT_EQ(result[3], 32.0);
    ASSERT_EQ(result[4], 86.0);
  }
}


TEST(CSRMatrix, ReadFromDisk)
{
  const std::string row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";

  auto spmat = readCSRMatrix(row_file, col_file, val_file);

  const std::string          rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";
  const std::vector<real_wp> rhs_vec  = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));
  ManagedArray<real_wp> temp_array({}, "test_temp", rhs_vec.size());


  auto spmv_timer = Logger::get().timer("multiplying CSR matrix");
  spmat->matVecMultiply(rhs_array, temp_array);
}

TEST(LDUCSRMatrix, CreateFromCSR)
{
  auto spmat  = construct5x5TestMatrix();
  auto ldumat = std::make_unique<LDUCSRMatrix>(*spmat);

  ManagedArray<real_wp> ldu_prod({}, "ldu_prod", ldumat->cols(), 1.0);
  ManagedArray<real_wp> csr_prod({}, "csr_prod", spmat->cols(), 1.0);

  const int             nrows = spmat->rows();
  ManagedArray<real_wp> x({}, "x", nrows);

  auto x_arr = write_access_host(x);
  maDGForAllHost(irow, 0, nrows, { x_arr[irow] = irow + 1; });
  // set x={1,2,3,4,5}

  ldumat->matVecMultiply(x, ldu_prod);
  spmat->matVecMultiply(x, csr_prod);

  auto csr_result = read_access_host(csr_prod);
  ASSERT_EQ(csr_result[0], 49.0);
  ASSERT_EQ(csr_result[1], 18.0);
  ASSERT_EQ(csr_result[2], 56.0);
  ASSERT_EQ(csr_result[3], 32.0);
  ASSERT_EQ(csr_result[4], 86.0);

  auto ldu_result = read_access_host(ldu_prod);
  ASSERT_EQ(ldu_result[0], 49.0);
  ASSERT_EQ(ldu_result[1], 18.0);
  ASSERT_EQ(ldu_result[2], 56.0);
  ASSERT_EQ(ldu_result[3], 32.0);
  ASSERT_EQ(ldu_result[4], 86.0);
}

TEST(LDUCSRMatrix, CreateFromCSRAndCompareSpMV)
{
  const std::string          row_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_row";
  const std::string          col_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_col";
  const std::string          val_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db.mtx.csr_val";
  const std::string          rhs_file = std::string(TESTDIR) + "/sparse_matrices/poisson3Db_b.mtx.csr_rhs";
  const std::vector<real_wp> rhs_vec  = readArray<real_wp>(rhs_file);
  ManagedArray<real_wp> rhs_array({}, "test_rhs", Span<real_wp>(const_cast<real_wp*>(rhs_vec.data()), rhs_vec.size()));

  auto spmat  = readCSRMatrix(row_file, col_file, val_file);
  auto ldumat = std::make_unique<LDUCSRMatrix>(*spmat);

  ManagedArray<real_wp> csr_prod({}, "csr_prod", rhs_vec.size());
  ManagedArray<real_wp> ldu_prod({}, "ldu_prod", rhs_vec.size());

  spmat->matVecMultiply(rhs_array, csr_prod);
  ldumat->matVecMultiply(rhs_array, ldu_prod);

  auto csr_vec = read_access_host(csr_prod);
  auto ldu_vec = read_access_host(ldu_prod);

  for (int i = 0; i < csr_vec.size(); ++i) {
    ASSERT_LT(abs(csr_vec[i] - ldu_vec[i]), 1.e-10)
        << " i: " << i << " csr_vec: " << csr_vec[i] << " ldu_vec: " << ldu_vec[i];
  }
}


