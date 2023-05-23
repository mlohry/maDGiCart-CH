#pragma once


#include <limits>

#include "managed_array_owner.hpp"
#include "memory_handle.hpp"
#include "memory_manager.hpp"
#include "span.hpp"
#include "typedefs.hpp"

#include "exec_includes.hpp"
#include "utils/noncopyable.hpp"
#include "utils/uniquely_named_object.hpp"


template <typename T>
class ManagedArray final : public UniquelyNamedObject<ManagedArray<T>>, private NonCopyable {
 public:
  enum class SyncStatus { InSync, HostModified, DeviceModified };

  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, int_t size)
      : handle_(MemoryManager<T>::get().allocate(owner.objectName() + "::" + array_name, size)),
        owning_(true),
        sync_status_(SyncStatus::HostModified)
  {
    auto arr = writeHost();
    maDGForAllHost(i, 0, size, { arr[i] = 0.0; });
  }

  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, int_t size, T const_val)
      : handle_(MemoryManager<T>::get().allocate(owner.objectName() + "::" + array_name, size)),
        owning_(true),
        sync_status_(SyncStatus::HostModified)
  {
    auto arr = writeHost();
    maDGForAllHost(i, 0, size, { arr[i] = const_val; });
  }

  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, Span<T> initial_value)
      : handle_(MemoryManager<T>::get().allocate(owner.objectName() + "::" + array_name, initial_value.size())),
        owning_(true),
        sync_status_(SyncStatus::HostModified)
  {
    auto arr = writeHost();
    maDGForAllHost(i, 0, initial_value.size(), { arr[i] = initial_value[i]; });
  }


  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, T* host_data, T* dev_data, int_t size)
      : handle_(
            MemoryManager<T>::get().allocateView(owner.objectName() + "::" + array_name, size, host_data, dev_data)),
        owning_(false),
        sync_status_(SyncStatus::InSync)
  {
    auto arr = writeHost();
    maDGForAllHost(i, 0, size, { arr[i] = 0.0; });
  }


  ~ManagedArray()
  {
    if (owning_) {
      MemoryManager<T>::get().deallocate(handle());
    }
    else {
      MemoryManager<T>::get().deallocateView(handle());
    }
  }


  auto writeHost()
  {
    sync_status_ = SyncStatus::HostModified;
    return Span<T>(hostData(), size());
  }

  auto readHost() const
  {
    if (sync_status_ == SyncStatus::DeviceModified) {
      MemoryManager<T>::get().copyDeviceToHost(handle());
      sync_status_ = SyncStatus::InSync;
    }
    return Span<T>(hostData(), size());
  }

  auto readWriteHost()
  {
    if (sync_status_ == SyncStatus::DeviceModified) {
      MemoryManager<T>::get().copyDeviceToHost(handle());
    }
    sync_status_ = SyncStatus::HostModified;
    return Span<T>(hostData(), size());
  }

  auto writeDevice()
  {
    sync_status_ = SyncStatus::DeviceModified;
    return Span<T>(deviceData(), size());
  }

  auto readDevice() const
  {
    if (sync_status_ == SyncStatus::HostModified) {
      MemoryManager<T>::get().copyHostToDevice(handle());
      sync_status_ = SyncStatus::InSync;
    }
    return Span<T>(deviceData(), size());
  }

  auto readWriteDevice()
  {
    if (sync_status_ == SyncStatus::HostModified) {
      MemoryManager<T>::get().copyHostToDevice(handle());
    }
    sync_status_ = SyncStatus::DeviceModified;
    return Span<T>(deviceData(), size());
  }


  int_t size() const { return MemoryManager<T>::get().size(handle()); }

  void swap(ManagedArray<T>& other)
  {
    const MemoryHandle tmp_handle = handle_;
    handle_                       = other.handle_;
    other.handle_                 = tmp_handle;
    std::swap(owning_, other.owning_);
    std::swap(sync_status_, other.sync_status_);
  }


 private:
  MemoryHandle       handle_;
  bool               owning_;
  mutable SyncStatus sync_status_;

  MemoryHandle handle() const { return handle_; }

  // should just return spans here.
  //  T* data() const { return const_cast<T*>(MemoryManager<T>::get().data(handle())); }
  T* deviceData() const { return const_cast<T*>(MemoryManager<T>::get().deviceData(handle())); }
  T* hostData() const { return const_cast<T*>(MemoryManager<T>::get().hostData(handle())); }
};


using IndexArray  = ManagedArray<int_t>;
using ScalarArray = ManagedArray<real_wp>;
