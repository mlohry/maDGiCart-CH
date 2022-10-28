#pragma once


#include <limits>

#include "managed_array_owner.hpp"
#include "memory_handle.hpp"
#include "memory_manager.hpp"
#include "span.hpp"
#include "typedefs.hpp"

#include "utils/noncopyable.hpp"
#include "utils/uniquely_named_object.hpp"
#include "exec_includes.hpp"


template <typename T>
class ManagedArray final : public UniquelyNamedObject<ManagedArray<T>>, private NonCopyable {
 public:

  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, int_t size)
      : handle_(MemoryManager<T>::get().allocate(owner.objectName() + "::" + array_name, size)),
        sync_status_(SyncStatus::InSync)
  {
    auto arr = writeHost();
//    maDGForAllHost(i, 0, size, { arr[i] = std::numeric_limits<T>::signaling_NaN(); });
    maDGForAllHost(i, 0, size, { arr[i] = 0.0; });
  }

  ManagedArray(const ManagedArrayOwner& owner, const std::string& array_name, Span<T> initial_value)
      : handle_(MemoryManager<T>::get().allocate(owner.objectName() + "::" + array_name, initial_value.size())),
        sync_status_(SyncStatus::InSync)
  {
    auto arr = writeHost();
    maDGForAllHost(i, 0, initial_value.size(), { arr[i] = initial_value[i]; });
  }


  ~ManagedArray() { MemoryManager<T>::get().deallocate(handle()); }


  auto
  writeHost()
  {
    sync_status_ = SyncStatus::HostModified;
    return Span<T>(hostData(), size());
  }

  auto
  readHost() const
  {
    if (sync_status_ == SyncStatus::DeviceModified) {
      MemoryManager<T>::get().copyDeviceToHost(handle());
      sync_status_ = SyncStatus::InSync;
    }
    return Span<T>(hostData(), size());
  }

  auto
  readWriteHost()
  {
    if (sync_status_ == SyncStatus::DeviceModified) {
      MemoryManager<T>::get().copyDeviceToHost(handle());
    }
    sync_status_ = SyncStatus::HostModified;
    return Span<T>(hostData(), size());
  }

  auto
  writeDevice()
  {
    sync_status_ = SyncStatus::DeviceModified;
    return Span<T>(deviceData(), size());
  }

  auto
  readDevice() const
  {
    if (sync_status_ == SyncStatus::HostModified) {
      MemoryManager<T>::get().copyHostToDevice(handle());
      sync_status_ = SyncStatus::InSync;
    }
    return Span<T>(deviceData(), size());
  }

  auto
  readWriteDevice()
  {
    if (sync_status_ == SyncStatus::HostModified) {
      MemoryManager<T>::get().copyHostToDevice(handle());
    }
    sync_status_ = SyncStatus::DeviceModified;
    return Span<T>(deviceData(), size());
  }


  int_t
  size() const
  {
    return MemoryManager<T>::get().size(handle());
  }

  void
  deallocate()
  {
    MemoryManager<T>::get().deallocate(handle());
  }

  MemoryHandle
  handle() const
  {
    return handle_;
  }

  enum class SyncStatus { InSync, HostModified, DeviceModified };

  SyncStatus
  syncStatus() const
  {
    return sync_status_;
  }

 private:
  MemoryHandle       handle_;
  mutable SyncStatus sync_status_;


  // should just return spans here.
  //  T* data() const { return const_cast<T*>(MemoryManager<T>::get().data(handle())); }
  T*
  deviceData() const
  {
    return const_cast<T*>(MemoryManager<T>::get().deviceData(handle()));
  }
  T*
  hostData() const
  {
    return const_cast<T*>(MemoryManager<T>::get().hostData(handle()));
  }
};


using IndexArray  = ManagedArray<int_t>;
using ScalarArray = ManagedArray<real_wp>;
