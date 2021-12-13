#pragma once

// plain for loop wrapped in a lambda
#define maDGForAllLambdaLoopHost(indexname, startindex, endindex, expression) \
  [=](){ for (int indexname = startindex; indexname < endindex; ++indexname){ expression; }}();


#ifdef MADG_USE_RAJA
#include "RAJA/RAJA.hpp"
#define MADG_HOST_DEVICE RAJA_HOST_DEVICE

#define ReduceSumRealHostSeq RAJA::ReduceSum<RAJA::seq_reduce, real_t>
#define ReduceMinRealHostSeq RAJA::ReduceMin<RAJA::seq_reduce, real_t>
#define ReduceMaxRealHostSeq RAJA::ReduceMax<RAJA::seq_reduce, real_t>
#endif


#define read_access_host(arg) arg.readHost()
#define write_access_host(arg) arg.writeHost()
#define read_write_access_host(arg) arg.readWriteHost()

#define maDGForAllHost(indexname, startindex, endindex, expression) \
  RAJA::forall<RAJA::simd_exec>(RAJA::RangeSegment(startindex, endindex), [=] (int indexname) { expression });


#ifdef MADG_USE_GPU
#ifdef __NVCC__
#include "cublas_v2.h"
#define MADG_CUDA_BLOCK_SIZE 256
using device_exec_policy = RAJA::cuda_exec<MADG_CUDA_BLOCK_SIZE>;

#define maDGForAllDevice(indexname, startindex, endindex, expression) \
  RAJA::forall<device_exec_policy>(RAJA::RangeSegment(startindex, endindex), [=] __device__ (int indexname) { expression });

#define maDGForAll maDGForAllDevice

#define read_access(arg) arg.readDevice()
#define write_access(arg) arg.writeDevice()
#define read_write_access(arg) arg.readWriteDevice()

#define ReduceSumReal RAJA::ReduceSum<RAJA::cuda_reduce_atomic, real_t>
#define ReduceMinReal RAJA::ReduceMin<RAJA::cuda_reduce_atomic, real_t>
#define ReduceMaxReal RAJA::ReduceMax<RAJA::cuda_reduce_atomic, real_t>

#endif // __NVCC__
#endif // MADG_USE_GPU


#ifdef MADG_USE_SERIAL
#define maDGForAll maDGForAllHost
#define read_access(arg) arg.readHost()
#define write_access(arg) arg.writeHost()
#define read_write_access(arg) arg.readWriteHost()
#define ReduceSumReal ReduceSumRealHostSeq
#define ReduceMinReal ReduceMinRealHostSeq
#define ReduceMaxReal ReduceMaxRealHostSeq
#endif

#ifdef MADG_USE_OPENMP
#include <omp.h>
#define maDGForAll(indexname, startindex, endindex, expression) \
  RAJA::forall<RAJA::omp_parallel_for_exec>(RAJA::RangeSegment(startindex, endindex), [=] (int indexname) { expression });
#define read_access(arg) arg.readHost()
#define write_access(arg) arg.writeHost()
#define read_write_access(arg) arg.readWriteHost()
#define ReduceSumReal RAJA::ReduceSum<RAJA::omp_reduce, real_t>
#define ReduceMinReal RAJA::ReduceMin<RAJA::omp_reduce, real_t>
#define ReduceMaxReal RAJA::ReduceMax<RAJA::omp_reduce, real_t>
#endif
