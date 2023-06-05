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
#ifdef MADG_USE_HIP
#include <hip/hip_runtime.h>

#define HIPCHECK(cmd) myHipCheck(cmd, __FILE__, __LINE__)
void
myHipCheck(hipError_t e, const char* file, int line)
{
  if (e != hipSuccess) {
    printf("Failed: HIP error %s:%d '%s'\n", __FILE__, __LINE__, hipGetErrorString(e));
    exit(EXIT_FAILURE);
  }
}
#endif

template <typename BaseType>
MemoryManager<BaseType>::~MemoryManager()
{
  if (!host_data_.empty()) {
    std::cout << "Warning: " << this->objectName() << " has " << host_data_.size()
              << " host arrays still allocated on destruction.\n";

#ifndef NDEBUG
    std::cout << "look at my stack trace of size " << malloc_stack_trace_.size() << "\n";
    for (const auto& trace : malloc_stack_trace_) {
      std::cout << trace.second << std::endl;
    }
#endif
  }

  if (!borrowed_host_data_.empty()) {
    std::cout << "Warning: " << this->objectName() << " has " << borrowed_host_data_.size()
              << " borrowed host arrays still allocated on destruction.\n";
  }
#ifdef MADG_USE_GPU
  if (!device_data_.empty()) {
    std::cout << "Warning: " << this->objectName() << " has " << device_data_.size()
              << " device arrays still allocated on destruction.\n";
  }
  if (!borrowed_device_data_.empty()) {
    std::cout << "Warning: " << this->objectName() << " has " << borrowed_device_data_.size()
              << " borrowed device arrays still allocated on destruction.\n";
  }
#endif
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
#ifdef MADG_USE_HIP
  BaseType* gpu_data;
  HIPCHECK(hipMalloc((void**)&gpu_data, size * sizeof(BaseType)));
  device_data_.insert(std::pair<MapKeyType, BaseType*>(curr_index_, gpu_data));
#endif

  return handle;
}

template <typename BaseType>
MemoryHandle
MemoryManager<BaseType>::allocateView(
    const std::string& unique_name,
    size_t             size,
    BaseType*          host_data,
    BaseType*          dev_data)
{
  curr_index_++;

  borrowed_host_data_.insert(std::pair<MapKeyType, BaseType*>(curr_index_, host_data));

  auto handle = MemoryHandle{curr_index_, curr_generation_, unique_name};
  handles_.insert(std::pair<MapKeyType, MemoryHandle>(curr_index_, handle));

  data_sizes_.insert(std::pair<MapKeyType, int_t>(curr_index_, size * sizeof(BaseType)));

#ifdef MADG_USE_GPU
  borrowed_device_data_.insert(std::pair<MapKeyType, BaseType*>(curr_index_, dev_data));
#endif

  return handle;
}

template <typename BaseType>
void
MemoryManager<BaseType>::deallocateView(MemoryHandle handle)
{
  if (host_data_.find(handle.index) != host_data_.end()) {
    delete[] host_data_.at(handle.index);
    host_data_.erase(handle.index);
  }
  else if (borrowed_host_data_.find(handle.index) != borrowed_host_data_.end()) {
    borrowed_host_data_.erase(handle.index);
  }

  data_sizes_.erase(handle.index);

#ifdef MADG_USE_GPU
  if (device_data_.find(handle.index) != device_data_.end()) {
#ifdef MADG_USE_CUDA
    cudaErrorCheck(cudaFree(device_data_.at(handle.index)));
#endif
#ifdef MADG_USE_HIP
    HIPCHECK(hipFree(device_data_.at(handle.index)));
#endif
    device_data_.erase(handle.index);
  }
  else if (borrowed_device_data_.find(handle.index) != borrowed_device_data_.end()) {
    borrowed_device_data_.erase(handle.index);
  }
#endif  // MADG_USE_GPU

#ifndef NDEBUG
  malloc_stack_trace_.erase(handle.index);
#endif
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
  delete[] host_data_.at(handle.index);
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


#ifdef MADG_USE_HIP
  if (device_data_.find(handle.index) == device_data_.end()) {
    std::cout << "Warning: tried to deallocate non-existing device MemoryHandle " << handle << std::endl;
    return;
  }
  HIPCHECK(hipFree(device_data_.at(handle.index)));
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
    try {
      return borrowed_host_data_.at(handle.index);
    }
    catch (const std::out_of_range& oor) {
      std::cerr << "hostData lookup failure: " << oor.what() << "\n"
                << "MemoryHandle: " << handle << "\n";
      std::cout << "MemoryManager contents:\n" << *this << std::endl;
      abort();
    }
  }
}


template <typename BaseType>
BaseType*
MemoryManager<BaseType>::deviceData(MemoryHandle handle)
{
#ifdef MADG_USE_GPU
  try {
    return device_data_.at(handle.index);
  }
  catch (const std::out_of_range& oor) {
    try {
      return borrowed_device_data_.at(handle.index);
    }
    catch (const std::out_of_range& oor) {
      std::cerr << "deviceData lookup failure: " << oor.what() << "\n"
                << "MemoryHandle: " << handle << "\n";
      std::cout << "MemoryManager contents:\n" << *this << std::endl;
      abort();
    }
  }
#else
  std::cout << "Called MemoryManager::deviceData without GPU support." << std::endl;
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
#endif
#ifdef MADG_USE_HIP
  HIPCHECK(hipMemcpy(deviceData(handle), hostData(handle), data_sizes_.at(handle.index), hipMemcpyHostToDevice));
  n_host_to_device_transfers_++;
#endif
#if !defined(MADG_USE_GPU)
  std::cout << "Called MemoryManager::copyHostToDevice without GPU support." << std::endl;
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
#endif
#ifdef MADG_USE_HIP
  const int_t size_in_bytes = data_sizes_.at(handle.index);
  HIPCHECK(hipMemcpy(hostData(handle), deviceData(handle), size_in_bytes, hipMemcpyDeviceToHost));
  n_device_to_host_transfers_++;

  if (log_memory_transfers_) {
    Logger::get().logDeviceMemoryTransfer(
        "hipMemcpyDeviceToHost: array name: " + handle.name +
        " size: " + std::to_string(double(size_in_bytes) / (1024.0 * 1024.0)) + "MB");
  }
#endif
#if !defined(MADG_USE_GPU)
  std::cout << "Called MemoryManager::copyDeviceToHost without GPU support." << std::endl;
#endif
}

// explicit instantiations to make sure they're compiled by the device compiler
template class MemoryManager<float>;
template class MemoryManager<double>;
template class MemoryManager<int>;
