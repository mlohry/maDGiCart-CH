#include "memory_manager.hpp"

#ifdef MADG_USE_CUDA
#include <cuda.h>
#include <cuda_runtime_api.h>

#define cudaErrorCheck(ans)               \
  {                                       \
    gpuAssert((ans), __FILE__, __LINE__); \
  }
inline void
gpuAssert(cudaError_t code, const char* file, int line, bool abort = true)
{
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort)
      exit(code);
  }
}
#endif

template <typename BaseType>
MemoryManager<BaseType>::~MemoryManager()
{
  if (!host_data_.empty()) {
    std::cout << "Warning: " << this->objectName() << " has " << host_data_.size()
              << " arrays still allocated on destruction.\n";

#ifndef NDEBUG
    std::cout << "look at my stack trace of size " << malloc_stack_trace_.size() << "\n";
    for (const auto& trace : malloc_stack_trace_) {
      std::cout << trace.second << std::endl;
    }
#endif
  }
  else {
  }
}


template <typename BaseType>
MemoryHandle
MemoryManager<BaseType>::allocate(const std::string& unique_name, size_t size)
{
  curr_index_++;

  // auto* data = static_cast<BaseType*>(malloc(size * sizeof(BaseType)));
  BaseType* data = new BaseType[size];
  host_data_.insert(std::pair<MapKeyType, BaseType*>(curr_index_, data));

  auto handle = MemoryHandle{curr_index_, curr_generation_, unique_name};
  handles_.insert(std::pair<MapKeyType, MemoryHandle>(curr_index_, handle));

  data_sizes_.insert(std::pair<MapKeyType, int_t>(curr_index_, size * sizeof(BaseType)));

#ifdef MADG_USE_CUDA
  BaseType* gpu_data;
  cudaErrorCheck(cudaMalloc((void**)&gpu_data, size * sizeof(BaseType)));

  device_data_.insert(std::pair<MapKeyType, BaseType*>(curr_index_, gpu_data));
#endif

  return handle;
}


template <typename BaseType>
void
MemoryManager<BaseType>::deallocate(MemoryHandle handle)
{
  if (host_data_.find(handle.index) == host_data_.end()) {
    std::cout << "Warning: tried to deallocate non-existing host MemoryHandle " << handle << std::endl;
    return;
  }

  // free(host_data_.at(handle.index));
  delete [] host_data_.at(handle.index);
  host_data_.erase(handle.index);
  data_sizes_.erase(handle.index);

#ifdef MADG_USE_CUDA
  if (device_data_.find(handle.index) == device_data_.end()) {
    std::cout << "Warning: tried to deallocate non-existing device MemoryHandle " << handle << std::endl;
    return;
  }
  cudaErrorCheck(cudaFree(device_data_.at(handle.index)));
  device_data_.erase(handle.index);
#endif


#ifndef NDEBUG
  malloc_stack_trace_.erase(handle.index);
#endif
}



template <typename BaseType>
  BaseType*
  MemoryManager<BaseType>::hostData(MemoryHandle handle)
  {
    try {
      return host_data_.at(handle.index);
    }
    catch (const std::out_of_range& oor) {
      std::cerr << "lookup failure: " << oor.what() << "\n"
                << "MemoryHandle: " << handle << "\n";
      std::cout << "MemoryManager contents:\n" << *this << std::endl;
      abort();
    }
  }


template <typename BaseType>
  BaseType*
  MemoryManager<BaseType>::deviceData(MemoryHandle handle)
  {
#ifdef MADG_USE_CUDA
    try {
      return device_data_.at(handle.index);
    }
    catch (const std::out_of_range& oor) {
      std::cerr << "lookup failure: " << oor.what() << "\n"
                << "MemoryHandle: " << handle << "\n";
      std::cout << "MemoryManager contents:\n" << *this << std::endl;
      abort();
    }
#else
    std::cout << "Called MemoryManager::deviceData without CUDA support." << std::endl;
#endif
    return nullptr;
  }


template <typename BaseType>
  int_t
  MemoryManager<BaseType>::size(MemoryHandle handle) const
  {
    return data_sizes_.at(handle.index) / sizeof(BaseType);
  }

template <typename BaseType>
  void
  MemoryManager<BaseType>::copyHostToDevice(MemoryHandle handle)
  {
    profile();
#ifdef MADG_USE_CUDA
    cudaErrorCheck(
        cudaMemcpy(deviceData(handle), hostData(handle), data_sizes_.at(handle.index), cudaMemcpyHostToDevice));
    n_host_to_device_transfers_++;
#else
    std::cout << "Called MemoryManager::copyHostToDevice without CUDA support." << std::endl;
#endif
  }


template <typename BaseType>
  void
  MemoryManager<BaseType>::copyDeviceToHost(MemoryHandle handle)
  {
    profile();
#ifdef MADG_USE_CUDA
    const int_t size_in_bytes = data_sizes_.at(handle.index);
    cudaErrorCheck(cudaMemcpy(hostData(handle), deviceData(handle), size_in_bytes, cudaMemcpyDeviceToHost));
    n_device_to_host_transfers_++;


    if (log_memory_transfers_) {
      Logger::get().logDeviceMemoryTransfer(
          "cudaMemcpyDeviceToHost: array name: " + handle.name +
          " size: " + std::to_string(double(size_in_bytes) / (1024.0 * 1024.0)) + "MB");
    }

#else
    std::cout << "Called MemoryManager::copyDeviceToHost without CUDA support." << std::endl;
#endif
  }

// explicit instantiations to make sure they're compiled by the device compiler
template class MemoryManager<float>;
template class MemoryManager<double>;
template class MemoryManager<int>;
