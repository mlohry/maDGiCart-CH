#pragma once

#include "strided_array.hpp"
#include "array_set.hpp"
#include "typedefs.hpp"

template <class T>
class StridedArraySet : public ArraySet<StridedArray<T>>
{
 public:
  StridedArraySet(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      int_t                     nvecs,
      const StridedMatrixShape& shape)
      : ArraySet<StridedArray<T>>(owner, array_name, nvecs, shape)
  {
  }

  virtual ~StridedArraySet() = default;
};
