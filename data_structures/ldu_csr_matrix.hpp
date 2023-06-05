#pragma once

#include <memory>

#include "csr_matrix.hpp"
#include "sparse_matrix.hpp"

class LDUCSRMatrix : public ManagedArrayOwner, public SparseMatrix {
 public:
  LDUCSRMatrix(const CSRMatrix& full_csr_matrix);

  int    rows() const override { return nrows_; }
  int    cols() const override { return ncols_; }
  void   matVecMultiply(const ManagedArray<real_wp>& x, ManagedArray<real_wp>& b) const override;
  double rmsNorm(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const override;

  void jacobiSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const override;
  void jacobiSmoothCacheDiagonalInverse(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const;
  void jacobiSmoothBLASCacheDiagonalInverse(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight)
      const;

  void forwardGaussSeidelSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const;
  void backwardGaussSeidelSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const;
  void symmetricGaussSeidelSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const;

  double rmsNormDiagScaledRHS(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const;

  void printRow(int irow) const;

  void updateDiags() const;

 private:
  const int nrows_;
  const int ncols_;

  std::unique_ptr<CSRMatrix>             lower_;
  std::unique_ptr<CSRMatrix>             upper_;
  std::unique_ptr<ManagedArray<real_wp>> diag_;
  std::unique_ptr<ManagedArray<real_wp>> diag_inv_;

  std::unique_ptr<ManagedArray<real_wp>> tmp_work_vector_;
  std::unique_ptr<ManagedArray<real_wp>> tmp_work_vector2_;

  mutable bool initialized_for_unit_diagonal_ = false;
};
