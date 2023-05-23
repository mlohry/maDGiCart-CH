#pragma once

#include <exception>

#include "memory_handle.hpp"
#include "memory_manager.hpp"
#include "strided_array.hpp"
#include "typedefs.hpp"

#include "utils/uniquely_named_object.hpp"


template <typename T, typename MapType>
class ReferenceMatrixAccessor {
 public:
  ReferenceMatrixAccessor(T* const data, int_t rows, int_t cols)
      : data_(data), rows_(rows), cols_(cols), map_(data_, rows_, cols_)
  {
  }

  MADG_HOST_DEVICE inline const MapType& operator[](int_t /*ignored*/) const { return map_; }
  MADG_HOST_DEVICE inline MapType&       operator[](int_t /*ignored*/) { return map_; }

  MADG_HOST_DEVICE inline const T* data() const { return data_; }

  MADG_HOST_DEVICE inline auto const rows() const { return rows_; }
  MADG_HOST_DEVICE inline auto const cols() const { return cols_; }

 private:
  T* const    data_;
  const int_t rows_;
  const int_t cols_;
  MapType     map_;  // redundant after making this a constructed value
};


class ReferenceMatrix : public UniquelyNamedObject<ReferenceMatrix> {
 public:
  using T = real_wp;

  ReferenceMatrix(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      const StridedMatrixShape& shape)
      : array_(owner, array_name, shape.rows_per_cell * shape.cols_per_cell),
        rows_(shape.rows_per_cell),
        cols_(shape.cols_per_cell)
  {
  }

  ReferenceMatrix(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      const StridedMatrixShape& shape,
      Span<T> initial_value)
      : array_(owner, array_name, initial_value),
        rows_(shape.rows_per_cell),
        cols_(shape.cols_per_cell)
  {
  }



 	//  ReferenceMatrix(int_t rows, int_t cols)
  //      : array_(this->objectName(), "UnnamedReferenceMatrix", rows * cols), rows_(rows),
  //      cols_(cols)
  //  {
  //  }

  //  ReferenceMatrix() : ReferenceMatrix(0, 0) {}

  virtual ~ReferenceMatrix() {}

  using MutableMapType =
      Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;
  using ConstMapType =
      Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>;


  auto readHost() const
  {
    return ReferenceMatrixAccessor<T, ConstMapType>(asArray().readHost().data(), rows_, cols_);
  }
  auto writeHost()
  {
    return ReferenceMatrixAccessor<T, MutableMapType>(asArray().writeHost().data(), rows_, cols_);
  }
  auto readWriteHost()
  {
    return ReferenceMatrixAccessor<T, MutableMapType>(
        asArray().readWriteHost().data(), rows_, cols_);
  }

  auto readDevice() const
  {
    return ReferenceMatrixAccessor<T, ConstMapType>(asArray().readDevice().data(), rows_, cols_);
  }
  auto writeDevice()
  {
    return ReferenceMatrixAccessor<T, MutableMapType>(asArray().writeDevice().data(), rows_, cols_);
  }
  auto readWriteDevice()
  {
    return ReferenceMatrixAccessor<T, MutableMapType>(
        asArray().readWriteDevice().data(), rows_, cols_);
  }


  const ManagedArray<T>& asArray() const { return array_; }
  ManagedArray<T>&       asArray() { return array_; }

  int_t rows() const { return rows_; }
  int_t cols() const { return cols_; }

 private:
  ManagedArray<T> array_;

  int_t rows_;
  int_t cols_;
};


class ReferenceArray : public ReferenceMatrix
{
 public:
  ReferenceArray(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      const StridedMatrixShape& shape)
      : ReferenceMatrix(owner, array_name, shape)
  {
    assert(shape.cols_per_cell==1);
  }
};
