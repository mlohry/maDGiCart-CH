#pragma once

#include "managed_array.hpp"
#include "managed_array_owner.hpp"

class SparseMatrix {
 public:
  virtual int    rows() const                                                                                = 0;
  virtual int    cols() const                                                                                = 0;
  virtual void   matVecMultiply(const ManagedArray<real_wp>& x, ManagedArray<real_wp>& b) const              = 0;
  virtual double rmsNorm(const ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b) const               = 0;
  virtual void   jacobiSmooth(ManagedArray<real_wp>& x, const ManagedArray<real_wp>& b, double weight) const = 0;
};