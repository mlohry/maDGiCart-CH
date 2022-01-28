#pragma once

#include "managed_array.hpp"
#include "utils/noncopyable.hpp"
#include "utils/uniquely_named_object.hpp"


MADG_HOST_DEVICE inline int
get1Dindex(int i, int j, int nhalo, int njhalo)
{
  return (i + nhalo) * njhalo + nhalo + j;
}


template <typename T>
class ManagedArray2DAccessor {
 public:
  ManagedArray2DAccessor(Span<T> span1d, int nhalo_layers, int nihalo, int njhalo)
      : span1d_(span1d), nhalo_layers_(nhalo_layers), nihalo_(nihalo), njhalo_(njhalo)
  {
  }

  MADG_HOST_DEVICE inline T& operator()(int i, int j) const { return span1d_[idx1d(i, j)]; }

  MADG_HOST_DEVICE inline void getIJ(int idx, int& i, int& j) const
  {
    i = idx / njhalo_ - nhalo_layers_;
    j = idx % njhalo_ - nhalo_layers_;
  }

 private:
  Span<T> span1d_;
  const int nhalo_layers_;
  const int nihalo_;
  const int njhalo_;

  MADG_HOST_DEVICE inline int idx1d(int i, int j) const { return get1Dindex(i, j, nhalo_layers_, njhalo_); }
};


template <typename T>
class ManagedArray2D : private NonCopyable {
 public:
  ManagedArray2D(int ni, int nj, int nhalo)
      : ni_(ni), nj_(nj), nhalo_(nhalo), nihalo_(ni + 2 * nhalo), njhalo_(nj + 2 * nhalo)
  {
  }

  virtual ~ManagedArray2D() = default;


  auto readWriteHost() { return ManagedArray2DAccessor<T>(asArray().readWriteHost(), nhalo_, nihalo_, njhalo_); }

  auto writeHost() { return ManagedArray2DAccessor<T>(asArray().writeHost(), nhalo_, nihalo_, njhalo_); }

  auto readHost() const { return ManagedArray2DAccessor<T>(asArray().readHost(), nhalo_, nihalo_, njhalo_); }


  auto readWriteDevice()
  {
    return ManagedArray2DAccessor<T>(asArray().readWriteDevice(), nhalo_, nihalo_, njhalo_);
  }

  auto writeDevice() { return ManagedArray2DAccessor<T>(asArray().writeDevice(), nhalo_, nihalo_, njhalo_); }

  auto readDevice() const { return ManagedArray2DAccessor<T>(asArray().readDevice(), nhalo_, nihalo_, njhalo_); }


  virtual const ManagedArray<T>& asArray() const = 0;
  virtual ManagedArray<T>&       asArray()       = 0;

 private:
  const int ni_;
  const int nj_;
  const int nhalo_;
  const int nihalo_;
  const int njhalo_;
};


template <typename T>
class ManagedArray2DOwning final : public ManagedArray2D<T> {
 public:
  ManagedArray2DOwning(const ManagedArrayOwner& owner, const std::string& array_name, int_t ni, int_t nj, int_t nhalo)
      : ManagedArray2D<T>(ni, nj, nhalo),
        array_(owner, array_name, (ni+2*nhalo)*(nj+2*nhalo))
  {
  }

  ~ManagedArray2DOwning() = default;


  const ManagedArray<T>& asArray() const override { return array_; }
  ManagedArray<T>&       asArray() override { return array_; }


 private:
  ManagedArray<T> array_;
};


template <typename T>
class ManagedArray2DNonOwning final : public ManagedArray2D<T> {
 public:
  ManagedArray2DNonOwning(ManagedArray<T>& src_array, int_t ni, int_t nj, int_t nhalo)
      :
      ManagedArray2D<T>(ni, nj, nhalo),
      borrowed_array_(src_array)
  {
  }

  ~ManagedArray2DNonOwning() = default;

  const ManagedArray<T>& asArray() const override { return borrowed_array_; }
  ManagedArray<T>&       asArray() override { return borrowed_array_; }


 private:
  ManagedArray<T>& borrowed_array_;
};
