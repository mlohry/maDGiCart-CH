#include "machine.hpp"

#include <boost/algorithm/string.hpp>
#include <fstream>
#include <vector>

#ifdef MADG_USE_OPENMP
#include <omp.h>
#endif

#ifdef MADG_USE_CUDA
#include <cuda_runtime.h>
#endif
#ifdef MADG_USE_HIP
#include <hip/hip_runtime.h>
#endif

#include <iostream>


Machine::Machine()
{
#ifdef __NVCC__
  // Set the device for this process in a multi-gpu environment.
  // this is a naive round robin
  //  const int ndevices = getDevicesPerNode();
  //  const int deviceid = Parallel::Get().rank() % ndevices;
  //  cudaSetDevice(deviceid);
#endif
}


std::string
Machine::getProcessorModel() const
{
  std::string cpu_name = "not found";

  std::string   line;
  std::ifstream cpuinfo("/proc/cpuinfo");
  if (cpuinfo.is_open()) {
    while (std::getline(cpuinfo, line)) {
      if (line.substr(0, 10) == "model name") {
        std::vector<std::string> words;
        boost::split(words, line, [](char c) { return c == ' '; });
        if (words.size() > 2) {
          cpu_name = words[2];
          for (size_t i = 3; i < words.size(); ++i) {
            cpu_name += " ";
            cpu_name += words[i];
          }
        }
        break;
      }
    }
  }

  return cpu_name;
}


std::string
Machine::getHostname() const
{
  char host_name[HOST_NAME_MAX];
  gethostname(host_name, HOST_NAME_MAX);  // posix system function, can return null
  return std::string(host_name);
}


std::string
Machine::getUsername() const
{
  char* user_name = std::getenv("USER");  // posix system function, can return null
  if (user_name != nullptr) {
    return std::string(user_name);
  }
  return "unknown";
}


std::string
Machine::getDeviceModel() const
{
#ifdef MADG_USE_CUDA
  int dev = 0;
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, dev);
  return std::string(deviceProp.name);
#endif
#ifdef MADG_USE_HIP
  int dev = 0;
  hipDeviceProp_t deviceProp;
  hipGetDeviceProperties(&deviceProp, dev);
  return std::string(deviceProp.name);
#endif
  return "none";
}


int
Machine::getDevicesPerNode() const
{
  int deviceCount = 0;
#ifdef MADG_USE_CUDA
  cudaGetDeviceCount(&deviceCount);
#endif
#ifdef MADG_USE_HIP
  hipGetDeviceCount(&deviceCount);
#endif
  return deviceCount;
}


std::string
Machine::getOpenMPInfo() const
{
#ifdef MADG_USE_OPENMP
  const int nthreads = omp_get_max_threads();
  return std::to_string(nthreads);
#endif
  return "Not compiled with OpenMP support";
}
