cmake_minimum_required(VERSION 3.2.2)
project(maDGiCart)

include(ExternalProject)

find_package(Boost REQUIRED COMPONENTS regex program_options log log_setup filesystem thread)

set(CXX_COMPILER_STRING "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
add_definitions(-DCXX_COMPILER_STRING="${CXX_COMPILER_STRING}")
add_definitions(-DBOOST_LOG_DYN_LINK)
add_definitions(-DTESTDIR=\"${CMAKE_SOURCE_DIR}/testing/\")
add_definitions(-DSOURCEDIR=\"${CMAKE_SOURCE_DIR}\")
add_definitions(-DBOOST_MATH_DISABLE_FLOAT128)
add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS) # remove annoying boost pragma compiler message


set(MADG_USE_RAJA ON)
add_definitions(-DMADG_USE_RAJA)

option(MADG_USE_CUDA "Build with CUDA GPU support" OFF)
option(MADG_USE_HIP "Build with HIP GPU support" OFF)
option(MADG_USE_GPU "Build with GPU support" OFF)
option(MADG_USE_OPENMP "Build with OpenMP support" OFF)
option(MADG_USE_SERIAL "Build with serial support" ON)
option(MADG_USE_SINGLE_PRECISION "Build with single precision instead of double" OFF)

if(MADG_USE_SINGLE_PRECISION)
    add_definitions(-DMADG_USE_SINGLE_PRECISION)
endif()

if (MADG_USE_CUDA)
    set(MADG_USE_GPU ON)
    add_definitions(-DMADG_USE_CUDA)
    enable_language(CUDA)
    find_package(CUDA REQUIRED)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcudafe --diag_suppress=177  --Wext-lambda-captures-this -arch=sm_35 --expt-extended-lambda --expt-relaxed-constexpr -Wno-deprecated-gpu-targets -Xcompiler -Wno-non-template-friend,-Wall")
endif(MADG_USE_CUDA)

if (MADG_USE_HIP)
    set(MADG_USE_GPU ON)
    add_definitions(-DMADG_USE_HIP)
    list (APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)
    find_package(hip REQUIRED)
    find_package(rocsparse REQUIRED)
    add_definitions(-DMADG_USE_ROCSPARSE)

endif(MADG_USE_HIP)

if (MADG_USE_GPU)
  if (MADG_USE_OPENMP)
    message(FATAL_ERROR "Only one of MADG_USE_GPU or MADG_USE_OPENMP can be enabled")
  endif()
endif()

if (MADG_USE_GPU)
  add_definitions(-DMADG_USE_GPU)
  set(MADG_USE_SERIAL OFF)
  MESSAGE(STATUS "MADG_USE_GPU specified, setting MADG_USE_SERIAL Off.")
endif()

if (MADG_USE_OPENMP)
  add_definitions(-DMADG_USE_OPENMP)
  set(MADG_USE_SERIAL OFF)
  MESSAGE(STATUS "MADG_USE_OPENMP specified, setting MADG_USE_SERIAL Off.")
  find_package(OpenMP REQUIRED)
endif()

if(MADG_USE_SERIAL)
  add_definitions(-DMADG_USE_SERIAL)
  MESSAGE(STATUS "MADG_USE_SERIAL is On.")
endif()


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if (CMAKE_BUILD_TYPE STREQUAL "Release")
    add_definitions(-DNDEBUG)
endif()


set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror -Wno-error=unused-function -Wno-error=sign-compare -Wno-int-in-bool-context -Werror=reorder -Wno-unused-result -fno-omit-frame-pointer ${STD_FLAG} ${FP_PRECISION_FLAG} ${OpenMP_CXX_FLAGS}")

#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -Wall -Werror -Wno-error=maybe-uninitialized ${STD_FLAG} ${FP_PRECISION_FLAG} -march=native ${OpenMP_CXX_FLAGS} -ftree-vectorize -ftree-vectorizer-verbose=1")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g ${STD_FLAG} ${FP_PRECISION_FLAG} -march=native ${OpenMP_CXX_FLAGS}")


set(EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external_install)
add_subdirectory(external)




include_directories(
        ${CMAKE_SOURCE_DIR}
        ${GTEST_INCLUDE_DIR}
        ${EIGEN3_INCLUDE_DIR}
        ${Boost_INCLUDE_DIRS}
        ${TABULATE_INCLUDE_DIR}
        ${CUDA_INCLUDE_DIRS}
        ${CUDA_MATH_INCLUDE_DIR}
        ${RAJA_INCLUDE_DIRS}
        ${EXTERNAL_INSTALL_DIR}/include
        ${TINYXML2_INCLUDE_DIR}
        ${PETSC_INCLUDES}
        ${ROC_LIB_INCLUDES}
)

link_directories(${EXTERNAL_INSTALL_DIR}/lib)

set(LIBTYPE "SHARED")

add_subdirectory(file_io)
add_subdirectory(logger)
add_subdirectory(data_structures)
add_subdirectory(governing_equations)
add_subdirectory(initialization)
add_subdirectory(parallel)
add_subdirectory(program_options)
add_subdirectory(spatial_discretization)
add_subdirectory(time_stepping)


# Create compile_commands.json in the build directory.
# Useful for clang-tidy, VSCode's cquery, and other analysis tools.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if (${LIBTYPE} STREQUAL "SHARED")
    set(LINKPREFIX "-Wl,--no-as-needed")
    set(LINKSUFFIX "-Wl,--as-needed")
else()
    set(LINKPREFIX "-Wl,--whole-archive")
    set(LINKSUFFIX "-Wl,--no-whole-archive")
endif()

set( maDG_library_links
        data_structures
        ${LINKPREFIX} governing_equations ${LINKSUFFIX}
        parallel
        discretization
        initialization
        logger
        ${LINKPREFIX} time_stepping ${LINKSUFFIX}
        ${LINKPREFIX} program_options ${LINKSUFFIX}
        )

add_executable(maDGiCart
        maDGiCart_main.cpp
        typedefs.hpp)

target_link_libraries(maDGiCart ${maDG_library_links} discretization RAJA)

enable_testing()
add_subdirectory(testing)