#pragma once

#include "exec_includes.hpp"
#include "managed_array.hpp"
#include "utils/noncopyable.hpp"
#include "utils/uniquely_named_object.hpp"

MADG_HOST_DEVICE inline int
get1DindexFrom3D(int i, int j, int k, int nhalo, int njhalo, int nkhalo)
{
  return (i + nhalo) * (njhalo * nkhalo) + (j + nhalo) * nkhalo + (k + nhalo);
}


template <typename T>
class ManagedArray3DAccessor {
 public:
  ManagedArray3DAccessor(T* const data, int nhalo_layers, int nihalo, int njhalo, int nkhalo)
      : data_(data), nhalo_layers_(nhalo_layers), nihalo_(nihalo), njhalo_(njhalo), nkhalo_(nkhalo)
  {
  }

#ifndef NDEBUG
  MADG_HOST_DEVICE inline T& operator()(int i, int j, int k) const {
    assert(i>= -nhalo_layers_);
    assert(i < nihalo_);
    assert(j>= -nhalo_layers_);
    assert(j < njhalo_);
    assert(k>= -nhalo_layers_);
    assert(k < nkhalo_);

    return data_[idx1d(i, j, k)];
  }
#else
  MADG_HOST_DEVICE inline T& operator()(int i, int j, int k) const { return data_[idx1d(i, j, k)]; }
#endif

  MADG_HOST_DEVICE inline void getIJK(int idx, int& i, int& j, int& k) const
  {
    i = idx / (njhalo_ * nkhalo_) - nhalo_layers_;
    j = (idx / (nkhalo_)) % njhalo_ - nhalo_layers_;
    k = idx % ( nkhalo_) - nhalo_layers_;
  }

 private:
  T* const  data_;
  const int nhalo_layers_;
  const int nihalo_;
  const int njhalo_;
  const int nkhalo_;

  MADG_HOST_DEVICE inline int idx1d(int i, int j, int k) const
  {
    return get1DindexFrom3D(i, j, k, nhalo_layers_, njhalo_, nkhalo_);
  }
};


template <typename T>
class ManagedArray3D : private NonCopyable {
 public:
  ManagedArray3D(int ni, int nj, int nk, int nhalo)
      : ni_(ni),
        nj_(nj),
        nk_(nk),
        nhalo_(nhalo),
        nihalo_(ni + 2 * nhalo),
        njhalo_(nj + 2 * nhalo),
        nkhalo_(nk + 2 * nhalo)
  {
  }

  virtual ~ManagedArray3D() = default;


  auto readWriteHost()
  {
    return ManagedArray3DAccessor<T>(asArray().readWriteHost().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }

  auto writeHost()
  {
    return ManagedArray3DAccessor<T>(asArray().writeHost().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }

  auto readHost() const
  {
    return ManagedArray3DAccessor<T>(asArray().readHost().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }


  auto readWriteDevice()
  {
    return ManagedArray3DAccessor<T>(asArray().readWriteDevice().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }

  auto writeDevice()
  {
    return ManagedArray3DAccessor<T>(asArray().writeDevice().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }

  auto readDevice() const
  {
    return ManagedArray3DAccessor<T>(asArray().readDevice().data(), nhalo_, nihalo_, njhalo_, nkhalo_);
  }


  virtual const ManagedArray<T>& asArray() const = 0;
  virtual ManagedArray<T>&       asArray()       = 0;

 private:
  const int ni_;
  const int nj_;
  const int nk_;
  const int nhalo_;
  const int nihalo_;
  const int njhalo_;
  const int nkhalo_;
};

template <typename T>
class ManagedArray3DOwning final : public ManagedArray3D<T> {
 public:
  ManagedArray3DOwning(
      const ManagedArrayOwner& owner,
      const std::string&       array_name,
      int_t                    ni,
      int_t                    nj,
      int_t                    nk,
      int_t                    nhalo)
      : ManagedArray3D<T>(ni, nj, nk, nhalo),
        array_(owner, array_name, (ni + 2 * nhalo) * (nj + 2 * nhalo) * (nk + 2 * nhalo))
  {
  }

  ~ManagedArray3DOwning() = default;


  const ManagedArray<T>& asArray() const override { return array_; }
  ManagedArray<T>&       asArray() override { return array_; }


 private:
  ManagedArray<T> array_;
};


template <typename T>
class ManagedArray3DNonOwning final : public ManagedArray3D<T> {
 public:
  ManagedArray3DNonOwning(ManagedArray<T>& src_array, int_t ni, int_t nj, int_t nk, int_t nhalo)
      : ManagedArray3D<T>(ni, nj, nk, nhalo), borrowed_array_(src_array)
  {
  }

  ~ManagedArray3DNonOwning() = default;

  const ManagedArray<T>& asArray() const override { return borrowed_array_; }
  ManagedArray<T>&       asArray() override { return borrowed_array_; }


 private:
  ManagedArray<T>& borrowed_array_;
};
