#pragma once

#include <exception>
#include <iostream>
#include <unordered_map>
#include "memory_handle.hpp"
#include "typedefs.hpp"

#include "logger/logger.hpp"
#include "logger/profiler.hpp"
#include "utils/uniquely_named_object.hpp"


template <typename BaseType>
class MemoryManager : public UniquelyNamedObject<MemoryManager<BaseType>> {

 public:
  using MapKeyType = int_t;

  static MemoryManager& get()
  {
    static MemoryManager instance;
    return instance;
  }

  // todo make allocate and deallocate private with friend class ManagedArray
  MemoryHandle allocate(const std::string& unique_name, size_t size);

  void deallocate(MemoryHandle handle);

  // todo these should return spans
  BaseType* hostData(MemoryHandle handle);

  BaseType* deviceData(MemoryHandle handle);

  int_t size(MemoryHandle handle) const;


  template <typename T>
  friend std::ostream& operator<<(std::ostream& os, const MemoryManager<T>& mempool);

  void copyHostToDevice(MemoryHandle handle);

  void copyDeviceToHost(MemoryHandle handle);


  bool log_memory_transfers_ = false;

 private:
  MemoryManager() {}

  ~MemoryManager();

  int_t curr_index_      = 0;
  int_t curr_generation_ = 0;

  unsigned long long n_host_to_device_transfers_ = 0;
  unsigned long long n_device_to_host_transfers_ = 0;

  std::unordered_map<MapKeyType, BaseType*>    host_data_;
  std::unordered_map<MapKeyType, int_t>        data_sizes_;
  std::unordered_map<MapKeyType, MemoryHandle> handles_;

#ifdef MADG_USE_CUDA
  std::unordered_map<MapKeyType, BaseType*> device_data_;
#endif

#ifndef NDEBUG
  std::unordered_map<MapKeyType, std::string> malloc_stack_trace_;
#endif
};


template <typename T>
std::ostream&
operator<<(std::ostream& os, const MemoryManager<T>& mempool)
{
  os << "MemoryManager size: " << mempool.host_data_.size() << "\n";
  for (const auto& kv : mempool.host_data_) {
    os << "Key: " << kv.first << "\nMemoryHandle: " << mempool.handles_.at(kv.first) << "\n";
  }
  return os;
}
