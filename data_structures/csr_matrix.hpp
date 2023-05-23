#pragma once

#include "managed_array.hpp"
#include "sparse_matrix.hpp"

#ifdef MADG_USE_ROCSPARSE
#include <rocsparse/rocsparse.h>
#endif

class CSRMatrix : public ManagedArrayOwner, public SparseMatrix {
 public:
  CSRMatrix(int_t rows, int_t cols, int_t nnonzeros);

  CSRMatrix(
      int                         nrow,
      int                         ncol,
      const std::vector<int>&     col,
      const std::vector<int>&     row,
      const std::vector<real_wp>& val);


  CSRMatrix(
      const CSRMatrix& source,
      real_wp*         host_values,
      int_t*           host_col_idx,
      int_t*           host_row_idx,
      real_wp*         dev_values,
      int_t*           dev_col_idx,
      int_t*           dev_row_idx);

  ~CSRMatrix();

  int    rows() const override { return nrows_; }
  int    cols() const override { return ncols_; }
  void   matVecMultiply(const ManagedArray<real_wp>& x, ManagedArray<real_wp>& b) const override;
  double rmsNorm(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const override;
  void   jacobiSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const override;

  void matVecMultiplyBLAS(const ManagedArray<real_wp>& x, ManagedArray<real_wp>& b, real_wp beta_in = 0) const;

  void jacobiSmoothCacheDiagonalInverse(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const;
  void jacobiSmoothBLASCacheDiagonalInverse(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight)
      const;


  const ManagedArray<real_wp>& values() const { return values_; }
  const ManagedArray<int_t>&   colIdx() const { return col_idx_; }
  const ManagedArray<int_t>&   rowIdx() const { return row_idx_; }

  ManagedArray<real_wp>& values()
  {
    diags_need_update_ = true;
    return const_cast<ManagedArray<real_wp>&>(static_cast<const CSRMatrix&>(*this).values());
  }
  ManagedArray<int_t>& colIdx()
  {
    diags_need_update_ = true;
    return const_cast<ManagedArray<int_t>&>(static_cast<const CSRMatrix&>(*this).colIdx());
  }
  ManagedArray<int_t>& rowIdx()
  {
    diags_need_update_ = true;
    return const_cast<ManagedArray<int_t>&>(static_cast<const CSRMatrix&>(*this).rowIdx());
  }

  enum class TriangularType { Upper, Lower };
  void scaleDiagonalAndInitTriangularSolve(const ManagedArray<real_wp>& diag_scale, TriangularType type);
  void scaleMatrixRows(const ManagedArray<real_wp>& diag_scale);
  void triangularSolve(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const;

  void printRow(int irow) const;

  real_wp diagonalDominance() const;

  void updateDiags() const;

 private:
  const int_t nrows_;
  const int_t ncols_;

  ManagedArray<real_wp>                          values_;
  ManagedArray<int_t>                            col_idx_;
  ManagedArray<int_t>                            row_idx_;
  mutable std::unique_ptr<ManagedArray<real_wp>> diag_inv_;
  mutable bool                                   diags_need_update_ = true;

  std::unique_ptr<ManagedArray<real_wp>> tmp_work_vector_;


#ifdef MADG_USE_ROCSPARSE
  rocsparse_handle    rocsparse_handle_;
  rocsparse_mat_descr spmv_descrA_;
  rocsparse_mat_info  spmv_mat_info_ = NULL;

  rocsparse_mat_descr                    csrsv_descrA_;
  rocsparse_mat_info                     csrsv_mat_info_ = NULL;

  std::unique_ptr<ManagedArray<real_wp>> csrsv_buffer_;
#endif

      void
       initCommon();
  void initSpMV();
};
