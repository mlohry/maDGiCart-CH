#pragma once

#include <vector>
#include <memory>
#include <cassert>

#include "managed_array_owner.hpp"


//template <template<class> class ArrayType, class ScalarType>
template <typename ArrayType>
class ArraySet
{
 public:
  ArraySet(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      int_t                     nvecs,
      int_t                     vecsize)
  {
    for (int_t i = 0; i < nvecs; ++i) {
      vecs_.push_back(
          std::make_unique<ArrayType>(owner, array_name + "Vec" + std::to_string(i), vecsize));
    }
  }


  virtual ~ArraySet() = default;

  using BaseType = ArrayType;


  const ArrayType& getVec(int_t i) const
  {
    assert(i >= 0);
    assert(i < nvecs());
    return *vecs_[i];
  }


  ArrayType& getVec(int_t i)
  {
    assert(i >= 0);
    assert(i < nvecs());
    return *vecs_[i];
  }


  int_t nvecs() const { return static_cast<int_t>(vecs_.size()); }


 private:
  std::vector<std::unique_ptr<ArrayType>> vecs_;
};
