#pragma once

#include "array_set.hpp"
#include "managed_array_owner.hpp"
#include "strided_array.hpp"


template <typename T>
class NodalScalarBase : public StridedArray<T> {
 public:
  NodalScalarBase(const ManagedArrayOwner& owner, const std::string& array_name, const StridedMatrixShape& shape)
      : StridedArray<T>(owner, array_name, shape)
  {
  }

  virtual ~NodalScalarBase() {}
};


class NodalScalarSet : public ArraySet<NodalScalarBase<real_wp>> {
 public:
  NodalScalarSet(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      int                       nvecs,
      const StridedMatrixShape& shape)
      : ArraySet<NodalScalarBase<real_wp>>(owner, array_name, nvecs, shape)
  {
  }
};
