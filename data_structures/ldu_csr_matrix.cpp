#include "ldu_csr_matrix.hpp"

#include <iostream>
#include <vector>


LDUCSRMatrix::LDUCSRMatrix(const CSRMatrix& full_csr_matrix)
    : nrows_(full_csr_matrix.rows()), ncols_(full_csr_matrix.cols())
{
  std::vector<real_wp> diags;
  std::vector<real_wp> upper_vals;
  std::vector<real_wp> lower_vals;
  std::vector<int>     upper_cols;
  std::vector<int>     lower_cols;
  std::vector<int>     upper_rows(nrows_ + 1, 0);
  std::vector<int>     lower_rows(nrows_ + 1, 0);

  auto col_idx = read_access_host(full_csr_matrix.colIdx());
  auto row_idx = read_access_host(full_csr_matrix.rowIdx());
  auto vals    = read_access_host(full_csr_matrix.values());

  for (int irow = 0; irow < nrows_; ++irow) {
    int nval_lower_in_row = 0;
    int nval_upper_in_row = 0;
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      const int col = col_idx[j];
      if (col == irow) {
        diags.push_back(vals[j]);
      }
      else if (col < irow) {
        lower_vals.push_back(vals[j]);
        lower_cols.push_back(col);
        nval_lower_in_row++;
      }
      else {
        upper_vals.push_back(vals[j]);
        upper_cols.push_back(col);
        nval_upper_in_row++;
      }
    }
    lower_rows[irow + 1] = lower_rows[irow] + nval_lower_in_row;
    upper_rows[irow + 1] = upper_rows[irow] + nval_upper_in_row;
  }


  diag_ = std::make_unique<ManagedArray<real_wp>>(
      *this, "ldu_csr_diag", Span<real_wp>(const_cast<real_wp*>(diags.data()), diags.size()));
  diag_inv_         = std::make_unique<ManagedArray<real_wp>>(*this, "ldu_csr_diag_inv", diags.size());
  lower_            = std::make_unique<CSRMatrix>(nrows_, ncols_, lower_cols, lower_rows, lower_vals);
  upper_            = std::make_unique<CSRMatrix>(nrows_, ncols_, upper_cols, upper_rows, upper_vals);
  tmp_work_vector_  = std::make_unique<ManagedArray<real_wp>>(*this, "ldu_tmp_work_vector", ncols_);
  tmp_work_vector2_ = std::make_unique<ManagedArray<real_wp>>(*this, "ldu_tmp_work_vector2", ncols_);

  updateDiags();
}

void
LDUCSRMatrix::updateDiags() const
{
  auto D        = read_access(*diag_);
  auto diag_inv = write_access(*diag_inv_);
  maDGForAll(i, 0, D.size(), {
    diag_inv[i] =
        1.0 /
        (std::abs(D[i]) > std::numeric_limits<real_wp>::epsilon() ? D[i] : std::numeric_limits<real_wp>::epsilon());
  });
}


void
LDUCSRMatrix::matVecMultiply(const ManagedArray<real_wp>& x_in, ManagedArray<real_wp>& b_out) const
{
  auto lower_col_idx = read_access(lower_->colIdx());
  auto lower_row_idx = read_access(lower_->rowIdx());
  auto lower_vals    = read_access(lower_->values());
  auto upper_col_idx = read_access(upper_->colIdx());
  auto upper_row_idx = read_access(upper_->rowIdx());
  auto upper_vals    = read_access(upper_->values());
  auto x             = read_access(x_in);
  auto b             = write_access(b_out);
  auto diag          = read_access(*diag_);

  maDGForAll(irow, 0, rows(), {
    real_wp row_sum = 0;
    for (int j = lower_row_idx[irow]; j < lower_row_idx[irow + 1]; ++j) {
      row_sum += lower_vals[j] * x[lower_col_idx[j]];
    }
    for (int j = upper_row_idx[irow]; j < upper_row_idx[irow + 1]; ++j) {
      row_sum += upper_vals[j] * x[upper_col_idx[j]];
    }
    row_sum += diag[irow] * x[irow];
    b[irow] = row_sum;
  });
}


double
LDUCSRMatrix::rmsNorm(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const
{
  matVecMultiply(x, *tmp_work_vector_);

  auto          lhs = read_access(*tmp_work_vector_);
  auto          rhs = read_access(b);
  ReduceSumReal squared_norm(0);

  maDGForAll(i, 0, lhs.size(), { squared_norm += pow(lhs[i] - rhs[i], 2.0); });

  return sqrt(squared_norm.get() / lhs.size());
}

double
LDUCSRMatrix::rmsNormDiagScaledRHS(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const
{
  matVecMultiply(x, *tmp_work_vector_);

  auto          lhs  = read_access(*tmp_work_vector_);
  auto          rhs  = read_access(b);
  auto          Dinv = read_access(*diag_inv_);
  ReduceSumReal squared_norm(0);

  maDGForAll(i, 0, lhs.size(), { squared_norm += pow(lhs[i] - Dinv[i] * rhs[i], 2.0); });

  return sqrt(squared_norm.get() / lhs.size());
}

void
LDUCSRMatrix::jacobiSmooth(ManagedArray<real_wp>& x_in_out, const ManagedArray<real_wp>& b_in, double weight) const
{
  lower_->matVecMultiply(x_in_out, *tmp_work_vector_);
  upper_->matVecMultiply(x_in_out, *tmp_work_vector2_);

  auto Lx = read_access(*tmp_work_vector_);
  auto Ux = read_access(*tmp_work_vector2_);

  auto          x = read_write_access(x_in_out);
  auto          D = read_access(*diag_);
  auto          b = read_access(b_in);
  const real_wp w = weight;

  if (w == 1.0) {
    // unweighted jacobi
    maDGForAll(i, 0, rows(), {
      const real_wp d =
          std::abs(D[i]) > std::numeric_limits<real_wp>::epsilon() ? D[i] : std::numeric_limits<real_wp>::epsilon();

      x[i] = real_wp(1) / d * (b[i] - Lx[i] - Ux[i]);
    });
  }
  else {
    // weighted jacobi
    maDGForAll(i, 0, rows(), {
      const real_wp d =
          std::abs(D[i]) > std::numeric_limits<real_wp>::epsilon() ? D[i] : std::numeric_limits<real_wp>::epsilon();

      x[i] = w * real_wp(1) / d * (b[i] - Lx[i] - Ux[i]) + (1.0 - w) * x[i];
    });
  }
}


void
LDUCSRMatrix::jacobiSmoothCacheDiagonalInverse(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  lower_->matVecMultiply(x_in_out, *tmp_work_vector_);
  upper_->matVecMultiply(x_in_out, *tmp_work_vector2_);

  auto Lx = read_access(*tmp_work_vector_);
  auto Ux = read_access(*tmp_work_vector2_);

  auto x     = read_write_access(x_in_out);
  auto b     = read_access(b_in);
  auto D_inv = read_access(*diag_inv_);

  const real_wp w = weight;

  if (w == 1.0) {
    // unweighted jacobi
    maDGForAll(i, 0, rows(), { x[i] = D_inv[i] * (b[i] - Lx[i] - Ux[i]); });
  }
  else {
    // weighted jacobi
    maDGForAll(i, 0, rows(), { x[i] = w * D_inv[i] * (b[i] - Lx[i] - Ux[i]) + (1.0 - w) * x[i]; });
  }
}

void
LDUCSRMatrix::jacobiSmoothBLASCacheDiagonalInverse(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  lower_->matVecMultiplyBLAS(x_in_out, *tmp_work_vector_);       // compute tmp = L*x
  upper_->matVecMultiplyBLAS(x_in_out, *tmp_work_vector_, 1.0);  // compute tmp = tmp + U*x = (L+U)*x

  auto LUx   = read_access(*tmp_work_vector_);
  auto b     = read_access(b_in);
  auto D_inv = read_access(*diag_inv_);
  auto x     = read_write_access(x_in_out);

  const real_wp w = weight;

  if (w == 1.0) {
    // unweighted jacobi
    maDGForAll(i, 0, rows(), {  // weighted jacobi
      x[i] = D_inv[i] * (b[i] - LUx[i]);
    });
  }
  else {
    maDGForAll(i, 0, rows(), {  // unweighted jacobi
      x[i] = w * D_inv[i] * (b[i] - LUx[i]) + (1.0 - w) * x[i];
    });
  }
}


void
LDUCSRMatrix::forwardGaussSeidelSmooth(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  if (!initialized_for_unit_diagonal_) {
    const_cast<LDUCSRMatrix*>(this)->lower_->scaleDiagonalAndInitTriangularSolve(
        *diag_inv_, CSRMatrix::TriangularType::Lower);
    const_cast<LDUCSRMatrix*>(this)->upper_->scaleDiagonalAndInitTriangularSolve(
        *diag_inv_, CSRMatrix::TriangularType::Upper);
    initialized_for_unit_diagonal_ = true;

    auto D = write_access(*diag_);
    maDGForAll(i, 0, D.size(), { D[i] = 1.0; });
  }

  // compute tmp = U*x
  upper_->matVecMultiplyBLAS(x_in_out, *tmp_work_vector_);

  {
    auto bmUx  = read_write_access(*tmp_work_vector_);
    auto b     = read_access(b_in);
    auto x     = read_access(x_in_out);
    auto D_inv = read_access(*diag_inv_);
    auto old_x = write_access(*tmp_work_vector2_);

    // compute (b-Ux) for lower triangular solve, and copy previous value to apply relaxation later
    maDGForAll(i, 0, bmUx.size(), {
      bmUx[i]  = D_inv[i] * b[i] - bmUx[i];
      old_x[i] = x[i];
    });
  }

  // solve x = L*^(-1)( b - Ux )
  lower_->triangularSolve(x_in_out, *tmp_work_vector_);

  if (weight != 1.0)  // apply relaxation
  {
    auto x     = read_write_access(x_in_out);
    auto old_x = read_access(*tmp_work_vector2_);
    maDGForAll(i, 0, x.size(), {  //
      x[i] = (1.0 - weight) * old_x[i] + weight * x[i];
    });
  }
}

void
LDUCSRMatrix::backwardGaussSeidelSmooth(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  if (!initialized_for_unit_diagonal_) {
    const_cast<LDUCSRMatrix*>(this)->lower_->scaleDiagonalAndInitTriangularSolve(
        *diag_inv_, CSRMatrix::TriangularType::Lower);
    const_cast<LDUCSRMatrix*>(this)->upper_->scaleDiagonalAndInitTriangularSolve(
        *diag_inv_, CSRMatrix::TriangularType::Upper);
    initialized_for_unit_diagonal_ = true;

    auto D = write_access(*diag_);
    maDGForAll(i, 0, D.size(), { D[i] = 1.0; });
  }

  // compute tmp = L*x
  lower_->matVecMultiplyBLAS(x_in_out, *tmp_work_vector_);

  {
    auto bmLx  = read_write_access(*tmp_work_vector_);
    auto b     = read_access(b_in);
    auto x     = read_access(x_in_out);
    auto D_inv = read_access(*diag_inv_);
    auto old_x = write_access(*tmp_work_vector2_);

    // compute (b-Lx) for upper triangular solve, and copy previous value to apply relaxation later
    maDGForAll(i, 0, bmLx.size(), {
      bmLx[i]  = D_inv[i] * b[i] - bmLx[i];
      old_x[i] = x[i];
    });
  }

  // solve x = U*^(-1)( b - Lx )
  upper_->triangularSolve(x_in_out, *tmp_work_vector_);

  if (weight != 1.0)  // apply relaxation
  {
    auto x     = read_write_access(x_in_out);
    auto old_x = read_access(*tmp_work_vector2_);
    maDGForAll(i, 0, x.size(), {  //
      x[i] = (1.0 - weight) * old_x[i] + weight * x[i];
    });
  }
}

void
LDUCSRMatrix::symmetricGaussSeidelSmooth(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  forwardGaussSeidelSmooth(x_in_out, b_in, weight);
  backwardGaussSeidelSmooth(x_in_out, b_in, weight);
}

void
LDUCSRMatrix::printRow(int irow) const
{
  auto lower_col_idx = read_access(lower_->colIdx());
  auto lower_row_idx = read_access(lower_->rowIdx());
  auto lower_vals    = read_access(lower_->values());
  auto upper_col_idx = read_access(upper_->colIdx());
  auto upper_row_idx = read_access(upper_->rowIdx());
  auto upper_vals    = read_access(upper_->values());
  auto diag          = read_access(*diag_);

  for (int j = lower_row_idx[irow]; j < lower_row_idx[irow + 1]; ++j) {
    std::cout << "(" << irow << "," << lower_col_idx[j] << "): " << lower_vals[j] << " ";
  }
  std::cout << "(" << irow << "," << irow << "): " << diag[irow] << " ";

  for (int j = upper_row_idx[irow]; j < upper_row_idx[irow + 1]; ++j) {
    std::cout << "(" << irow << "," << upper_col_idx[j] << "): " << upper_vals[j] << " ";
  }
  std::cout << "\n";
}