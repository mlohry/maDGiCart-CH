#include "csr_matrix.hpp"

#include "exec_includes.hpp"


#ifdef MADG_USE_HIP
#include <hip/hip_runtime_api.h>
#include <rocsparse/rocsparse.h>

#define HIP_CHECK(stat)                                                  \
  {                                                                      \
    if (stat != hipSuccess) {                                            \
      std::cerr << "Error: hip error in line " << __LINE__ << std::endl; \
      return -1;                                                         \
    }                                                                    \
  }

#define ROCSPARSE_CHECK(stat)                                                  \
  {                                                                            \
    if (stat != rocsparse_status_success) {                                    \
      std::cerr << "Error: rocsparse error in line " << __LINE__ << std::endl; \
      return -1;                                                               \
    }                                                                          \
  }

#endif

CSRMatrix::CSRMatrix(int_t rows, int_t cols, int_t nnonzeros)
    : nrows_(rows),
      ncols_(cols),
      values_(*this, "CSRValues", nnonzeros),
      col_idx_(*this, "CSRColIndex", nnonzeros),
      row_idx_(*this, "CSRRowIndex", rows + 1)
{
  initCommon();
}

CSRMatrix::~CSRMatrix()
{
#ifdef MADG_USE_HIP
  rocsparse_csrmv_clear(rocsparse_handle_, spmv_mat_info_);
#endif
}

void
CSRMatrix::initCommon()
{
  tmp_work_vector_ = std::make_unique<ManagedArray<real_wp>>(*this, "csr_tmp_work_vector", ncols_);
  initSpMV();
}


void
CSRMatrix::initSpMV()
{
#ifdef MADG_USE_HIP
  //  int ierr = ROCSPARSE_CHECK(rocsparse_create_handle(&rocsparse_handle_));
  rocsparse_create_handle(&rocsparse_handle_);
  hipDeviceProp_t devProp;
  int             device_id = 0;
  hipGetDevice(&device_id);
  hipGetDeviceProperties(&devProp, device_id);
  //  std::cout << "Device: " << devProp.name << std::endl;
  rocsparse_create_mat_descr(&spmv_descrA_);

  auto vals    = read_access(values());
  auto row_idx = read_access(rowIdx());
  auto col_idx = read_access(colIdx());

#ifdef MADG_USE_SINGLE_PRECISION
  rocsparse_scsrmv_analysis(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      cols(),
      values_.size(),
      spmv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      spmv_mat_info_);
#else
  rocsparse_dcsrmv_analysis(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      cols(),
      values_.size(),
      spmv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      spmv_mat_info_);
#endif
#endif
}


void
CSRMatrix::scaleDiagonalAndInitTriangularSolve(const ManagedArray<real_wp>& diag_scale, TriangularType type)
{
  scaleMatrixRows(diag_scale);
#ifdef MADG_USE_ROCSPARSE
  rocsparse_create_mat_descr(&csrsv_descrA_);
  if (type == TriangularType::Upper) {
    rocsparse_set_mat_fill_mode(csrsv_descrA_, rocsparse_fill_mode_upper);
  }
  else {
    rocsparse_set_mat_fill_mode(csrsv_descrA_, rocsparse_fill_mode_lower);
  }
  rocsparse_set_mat_diag_type(csrsv_descrA_, rocsparse_diag_type_unit);
  rocsparse_create_mat_info(&csrsv_mat_info_);
  auto vals    = read_access(values());
  auto row_idx = read_access(rowIdx());
  auto col_idx = read_access(colIdx());

  // Obtain required buffer size
  size_t buffer_size;
  rocsparse_dcsrsv_buffer_size(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      values_.size(),
      csrsv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      csrsv_mat_info_,
      &buffer_size);

  csrsv_buffer_ = std::make_unique<ManagedArray<real_wp>>(*this, "csrsv_buffer", buffer_size/sizeof(real_wp));

  // Analysis policy
  rocsparse_analysis_policy analysis_policy = rocsparse_analysis_policy_reuse;
  // Solve policy
  rocsparse_solve_policy solve_policy = rocsparse_solve_policy_auto;
  auto                   buffer       = write_access(*csrsv_buffer_);
  // Perform analysis step
  rocsparse_dcsrsv_analysis(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      values_.size(),
      csrsv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      csrsv_mat_info_,
      analysis_policy,
      solve_policy,
      buffer.data());

#else
  Logger::get().FatalMessage("CSRMatrix::initForTriangularSolve not implemented");
#endif
}

void CSRMatrix::triangularSolve(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const
{

#ifdef MADG_USE_ROCSPARSE
  const real_wp alpha = 1.0;
  auto vals    = read_access(values());
  auto row_idx = read_access(rowIdx());
  auto col_idx = read_access(colIdx());
  auto sol = read_write_access(x);
  auto rhs = read_access(b);
  auto                   buffer       = write_access(*csrsv_buffer_);
  rocsparse_solve_policy solve_policy = rocsparse_solve_policy_auto;
  rocsparse_dcsrsv_solve(rocsparse_handle_,
                         rocsparse_operation_none,
                         rows(),
                         values_.size(),
                         &alpha,
                         csrsv_descrA_,
                         vals.data(),
                         row_idx.data(),
                         col_idx.data(),
                         csrsv_mat_info_,
                         rhs.data(),
                         sol.data(),
                         solve_policy,
                         buffer.data());

#else
  Logger::get().FatalMessage("CSRMatrix::initForTriangularSolve not implemented");
#endif
}


namespace {

template <typename T>
void
copyArray(const ManagedArray<T>& in, ManagedArray<T>& out)
{
  auto vout = write_access(out);
  auto vin  = read_access(in);
  maDGForAll(i, 0, vin.size(), { vout[i] = vin[i]; });
}


template <typename T>
void
copyArray(const std::vector<T>& in, ManagedArray<T>& out)
{
  auto vout = write_access_host(out);
  maDGForAllHost(i, 0, in.size(), { vout[i] = in[i]; });
}

}  // namespace


CSRMatrix::CSRMatrix(
    int                         nrow,
    int                         ncol,
    const std::vector<int>&     col,
    const std::vector<int>&     row,
    const std::vector<real_wp>& val)
    : nrows_(nrow),
      ncols_(ncol),
      values_(*this, "CSRValues", val.size()),
      col_idx_(*this, "CSRColIndex", col.size()),
      row_idx_(*this, "CSRRowIndex", row.size())
{
  copyArray(col, col_idx_);
  copyArray(row, row_idx_);
  copyArray(val, values_);
  initCommon();
}


CSRMatrix::CSRMatrix(
    const CSRMatrix& source,
    real_wp*         host_values,
    int_t*           host_col_idx,
    int_t*           host_row_idx,
    real_wp*         dev_values,
    int_t*           dev_col_idx,
    int_t*           dev_row_idx)
    : nrows_(source.rows()),
      ncols_(source.cols()),
      values_(*this, "CSRValues", host_values, dev_values, source.values().size()),
      col_idx_(*this, "CSRColIndex", host_col_idx, dev_col_idx, source.colIdx().size()),
      row_idx_(*this, "CSRRowIndex", host_row_idx, dev_row_idx, source.rowIdx().size())
{
  copyArray(source.values(), values());
  copyArray(source.colIdx(), colIdx());
  copyArray(source.rowIdx(), rowIdx());
  initCommon();
}


void
CSRMatrix::matVecMultiply(const ManagedArray<real_wp>& x_in, ManagedArray<real_wp>& b_out) const
{
  auto col_idx = read_access(colIdx());
  auto row_idx = read_access(rowIdx());
  auto vals    = read_access(values());
  auto x       = read_access(x_in);
  auto b       = write_access(b_out);

  maDGForAll(irow, 0, rows(), {
    real_wp row_sum = 0;
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      row_sum += vals[j] * x[col_idx[j]];
    }
    b[irow] = row_sum;
  });
}

void
CSRMatrix::matVecMultiplyBLAS(const ManagedArray<real_wp>& x_in, ManagedArray<real_wp>& b_out, real_wp beta_in) const
{

#ifdef MADG_USE_ROCSPARSE
  auto col_idx = read_access(colIdx());
  auto row_idx = read_access(rowIdx());
  auto vals    = read_access(values());
  auto x       = read_access(x_in);
  auto b       = write_access(b_out);


#ifdef MADG_USE_SINGLE_PRECISION
  const float alpha = 1.0;
  const float beta  = beta_in;

  rocsparse_scsrmv(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      cols(),
      values_.size(),
      &alpha,
      spmv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      spmv_mat_info_,
      x.data(),
      &beta,
      b.data());
  hipDeviceSynchronize();
#else
  const double alpha = 1.0;
  const double beta  = beta_in;

  rocsparse_dcsrmv(
      rocsparse_handle_,
      rocsparse_operation_none,
      rows(),
      cols(),
      values_.size(),
      &alpha,
      spmv_descrA_,
      vals.data(),
      row_idx.data(),
      col_idx.data(),
      NULL,
      x.data(),
      &beta,
      b.data());
  hipDeviceSynchronize();
#endif


#else
  Logger::get().FatalMessage("CSRMatrix::matVecMultiplyBLAS only implemented for rocsparse.");
#endif
}

double
CSRMatrix::rmsNorm(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const
{
  matVecMultiply(x, *tmp_work_vector_);

  auto          lhs = read_access(*tmp_work_vector_);
  auto          rhs = read_access(b);
  ReduceSumReal squared_norm(0);

  maDGForAll(i, 0, lhs.size(), { squared_norm += pow(lhs[i] - rhs[i], 2.0); });

  return sqrt(squared_norm.get() / lhs.size());
}


void
CSRMatrix::jacobiSmooth(ManagedArray<real_wp>& x_in_out, const ManagedArray<real_wp>& b_in, double weight) const
{
  auto          col_idx = read_access(colIdx());
  auto          row_idx = read_access(rowIdx());
  auto          vals    = read_access(values());
  auto          x       = read_write_access(x_in_out);
  auto          b       = read_access(b_in);
  auto          x_tmp   = write_access(*tmp_work_vector_);
  const real_wp w       = weight;

  maDGForAll(irow, 0, rows(), {
    real_wp row_sum = 0;
    real_wp d;
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      const int col = col_idx[j];
      if (col != irow) {
        row_sum += vals[j] * x[col];
      }
      else {
        d = std::abs(vals[j]) > std::numeric_limits<real_wp>::epsilon() ? vals[j]
                                                                        : std::numeric_limits<real_wp>::epsilon();
      }
    }
    if (w == 1.0) {
      x_tmp[irow] = 1.0 / d * (b[irow] - row_sum);
    }
    else {
      x_tmp[irow] = (w / d) * (b[irow] - row_sum) + (1.0 - w) * x[irow];
    }
  });

  x_in_out.swap(*tmp_work_vector_);
}


void
CSRMatrix::jacobiSmoothCacheDiagonalInverse(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  if (diags_need_update_) {
    updateDiags();
    diags_need_update_ = false;
  }

  auto col_idx = read_access(colIdx());
  auto row_idx = read_access(rowIdx());
  auto vals    = read_access(values());
  auto x       = read_write_access(x_in_out);
  auto b       = read_access(b_in);
  auto x_tmp   = write_access(*tmp_work_vector_);
  auto d_inv   = read_access(*diag_inv_);

  const real_wp w = weight;

  maDGForAll(irow, 0, rows(), {
    real_wp row_sum = 0;
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      row_sum += vals[j] * x[col_idx[j]];
    }
    x_tmp[irow] = w * d_inv[irow] * b[irow] + (x[irow] - w * d_inv[irow] * row_sum);
  });

  x_in_out.swap(*tmp_work_vector_);
}


void
CSRMatrix::jacobiSmoothBLASCacheDiagonalInverse(
    ManagedArray<real_wp>&       x_in_out,
    const ManagedArray<real_wp>& b_in,
    double                       weight) const
{
  if (diags_need_update_) {
    updateDiags();
    diags_need_update_ = false;
  }

  matVecMultiplyBLAS(x_in_out, *tmp_work_vector_);

  auto x     = read_access(x_in_out);
  auto b     = read_access(b_in);
  auto x_tmp = read_write_access(*tmp_work_vector_);
  auto d_inv = read_access(*diag_inv_);

  const real_wp w = weight;

  maDGForAll(irow, 0, rows(), { x_tmp[irow] = w * d_inv[irow] * b[irow] + (x[irow] - w * d_inv[irow] * x_tmp[irow]); });

  x_in_out.swap(*tmp_work_vector_);
}


void
CSRMatrix::updateDiags() const
{
  diag_inv_     = std::make_unique<ManagedArray<real_wp>>(*this, "csr_diag_inv", rows());
  auto col_idx  = read_access(colIdx());
  auto row_idx  = read_access(rowIdx());
  auto vals     = read_access(values());
  auto diag_inv = write_access(*diag_inv_);

  maDGForAll(irow, 0, rows(), {
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      const int col = col_idx[j];
      if (col == irow) {
        diag_inv[irow] = 1.0 / (std::abs(vals[j]) > std::numeric_limits<real_wp>::epsilon()
                                    ? vals[j]
                                    : std::numeric_limits<real_wp>::epsilon());
      }
    }
  });
}


void
CSRMatrix::scaleMatrixRows(const ManagedArray<real_wp>& diag_scale)
{
  auto row_idx   = read_access(rowIdx());
  auto vals      = read_write_access(values());
  auto row_scale = read_access(diag_scale);

  maDGForAll(irow, 0, rows(), {
    for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
      vals[j] *= row_scale[irow];
    }
  });
}


real_wp
CSRMatrix::diagonalDominance() const
{
  ReduceMinReal diag_dominance(std::numeric_limits<real_wp>::max());

  auto col_idx = read_access(colIdx());
  auto row_idx = read_access(rowIdx());
  auto vals    = read_access(values());

  maDGForAll(irow, 0, row_idx.size() - 1, {
    const int row_start = row_idx[irow];
    const int row_end   = row_idx[irow + 1];

    real_wp diag_term     = 0;
    real_wp off_diag_term = 0;

    for (int val_idx = row_start; val_idx < row_end; ++val_idx) {
      const int icol = col_idx[val_idx];

      if (icol == irow) {
        diag_term = std::abs(vals[val_idx]);
      }
      else {
        off_diag_term += std::abs(vals[val_idx]);
      }
    }

    diag_dominance.min(diag_term / off_diag_term);
  });

  return diag_dominance.get();
}

void
CSRMatrix::printRow(int irow) const
{
  auto col_idx = read_access(colIdx());
  auto row_idx = read_access(rowIdx());
  auto vals    = read_access(values());


  for (int j = row_idx[irow]; j < row_idx[irow + 1]; ++j) {
    std::cout << "(" << irow << "," << col_idx[j] << "): " << vals[j] << " ";
  }

  std::cout << "\n";
}