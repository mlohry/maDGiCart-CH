#pragma once


#include "typedefs.hpp"
#include "exec_includes.hpp"


template <typename T>
class Span
{
 public:
  Span(T* data, int_t size) : data_(data), size_(size) {}

  MADG_HOST_DEVICE inline T& operator[](int_t i) const {
    #ifndef NDEBUG
    assert(i >= 0);
    assert(i < size_);
    #endif
    return data_[i];
  }

  int_t size() const { return size_; }
  T* data() const { return data_; }


  // iterator access
  T* begin() const { return data(); }
  T* end() const { return data() + size(); }

 private:
  T* data_;
  int_t size_;
};
