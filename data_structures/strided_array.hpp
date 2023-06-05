#pragma once

#include <exception>

#include "memory_handle.hpp"
#include "memory_manager.hpp"
#include "typedefs.hpp"

#include "utils/uniquely_named_object.hpp"

#include "managed_array.hpp"
#include "span.hpp"

#include <memory>

#include <Eigen/Dense>


/**
 *
 * @tparam T       Scalar type
 * @tparam MapType A non-owning Eigen::Map type of per-cell data providing GPU-friendly linear algebra semantics
 */
template <typename T, typename MapType>
class StridedArrayAccessor
{
 public:
  /**
   *
   * @param data   Pointer to start of data, either device or host
   * @param stride Usually number of DOFs per cell. Real code this would be variable to accommodate mixed grids
   * @param ncells
   */
  StridedArrayAccessor(T* const data, int_t stride, int_t ncells)
      : data_(data), stride_(stride), ncells_(ncells)
  {
  }

  /**
   *
   * @param icell Index of the cell
   * @return A column vector of data in the cell constructed from the computed offset
   */
  MADG_HOST_DEVICE inline MapType operator[](int_t icell) const
  {
    return MapType(&(data_[icell * stride_]), stride_);
  }

  /**
   *
   * @param icell  Index of the cell
   * @return Pointer to the start of the cell's data, in case one wants to write raw kernel code (e.g. for optimization)
   */
  MADG_HOST_DEVICE inline T* data(int_t icell) const
  {
    return &(data_[icell * stride_]);
  }

  MADG_HOST_DEVICE inline const T* data() const { return data_; }
  MADG_HOST_DEVICE inline T* data() { return data_; }

 private:
  T* const    data_;
  const int_t stride_;
  const int_t ncells_;
};

struct StridedMatrixShape {
  int_t rows_per_cell = -1;
  int_t cols_per_cell = -1;
  int_t ncells        = -1;
  int_t size() const { return rows_per_cell*cols_per_cell*ncells; }
};



// todo construct an equivalent non-owning class or logic for borrowing ManagedArray's.
template <class T>
class StridedArray : private NonCopyable
{
 public:
  using ScalarType = T;

  StridedArray(
      const ManagedArrayOwner&  owner,
      const std::string&        array_name,
      const StridedMatrixShape& shape)
      : owned_array_(
        std::make_unique<ManagedArray<T>>(
            owner,
            array_name,
            shape.rows_per_cell * shape.cols_per_cell * shape.ncells)),
        borrowed_array_(nullptr),
        stride_(shape.rows_per_cell * shape.cols_per_cell),
        ncells_(shape.ncells)
  {
  }

  // constructor for non-owning
  StridedArray(
      ManagedArray<T>& array,
      const StridedMatrixShape& shape)
      : owned_array_(nullptr),
        borrowed_array_(&array),
        stride_(shape.rows_per_cell * shape.cols_per_cell),
        ncells_(shape.ncells)
  {
  }
  

  virtual ~StridedArray() {}

  using MutableMapType = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>>;
  using ConstMapType   = Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>>;


  auto readWriteHost()
  {
    return StridedArrayAccessor<T, MutableMapType>(
        asArray().readWriteHost().data(), stride_, ncells_);
  }

  auto writeHost() {
    return StridedArrayAccessor<T, MutableMapType>(
        asArray().writeHost().data(), stride_, ncells_);
  }

  auto readHost() const {
    return StridedArrayAccessor<T, ConstMapType>(
        asArray().readHost().data(), stride_, ncells_);
  }

  auto readWriteDevice()
  {
    return StridedArrayAccessor<T, MutableMapType>(
        asArray().readWriteDevice().data(), stride_, ncells_);
  }

  auto writeDevice() {
    return StridedArrayAccessor<T, MutableMapType>(
        asArray().writeDevice().data(), stride_, ncells_);
  }

  auto readDevice() const {
    return StridedArrayAccessor<T, ConstMapType>(
        asArray().readDevice().data(), stride_, ncells_);
  }

  auto handle() const { return asArray().handle(); }

  int_t stride() const { return stride_; }
  int_t ncells() const { return ncells_; } // basically outerStride() / innerStride()
  int_t size() const { return stride() * ncells(); }

  /*
 * todo: only public functions should be (a) flat spans, (b) strided vector accessors
 * no data() access
 * */

  const ManagedArray<T>& asArray() const { if (owned_array_){ return *owned_array_; }
  return *borrowed_array_; }
  ManagedArray<T>& asArray() {  if (owned_array_){ return *owned_array_; }
  return *borrowed_array_; }

 private:
 
  std::unique_ptr<ManagedArray<T>> owned_array_;
  ManagedArray<T>* borrowed_array_; // todo move the notion of "borrowed" upstream to ManagedArray itself

  int_t stride_ = 0;
  int_t ncells_ = 0;
};

