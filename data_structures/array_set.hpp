
#pragma once

#include "strided_array.hpp"


//template <template<class> class ArrayType, class ScalarType>
template <typename ArrayType>
class ArraySet
{
 public:
  ArraySet(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      int_t                     nvecs,
      const StridedMatrixShape& shape)
  {
    for (int_t i = 0; i < nvecs; ++i) {
      vecs_.push_back(
          std::make_unique<ArrayType>(owner, array_name + "Vec" + std::to_string(i), shape));
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


  using ScalarType = typename ArrayType::ScalarType;

  int_t nvecs() const { return static_cast<int_t>(vecs_.size()); }


 private:
  std::vector<std::unique_ptr<ArrayType>> vecs_;
};
