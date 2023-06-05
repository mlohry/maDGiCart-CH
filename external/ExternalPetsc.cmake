#########################################################
message( STATUS "Building external Petsc project.")
#########################################################

if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(PETSC_OPT_FLAGS "-O3 -march=native")
    set(PETSC_DEBUGGING "no")
    set(PETSC_ARCH_FLAG "arch-linux2-c-opt")
else()
    set(PETSC_OPT_FLAGS "-O1 -g")
    set(PETSC_DEBUGGING "yes")
    set(PETSC_ARCH_FLAG "linux-gnu-c-debug")
endif()

if(MADG_USE_CUDA)
    set(PETSC_CUDAC_FLAG "--with-cudac=${CMAKE_CUDA_COMPILER}")
    set(PETSC_CUDA_FLAG "--with-cuda")
    set(PETSC_CUDA_ARCH_FLAG "--with-cuda-arch=35")
    set(PETSC_CUDA_DIR_FLAG "--with-cuda-dir=${CUDA_TOOLKIT_ROOT_DIR}")
else()
    set(PETSC_CUDAC_FLAG "")
    set(PETSC_CUDA_FLAG "")
    set(PETSC_CUDA_ARCH_FLAG "")
    set(PETSC_CUDA_DIR_FLAG "")
endif()

if (MADG_USE_SINGLE_PRECISION)
    set(PETSC_PRECISION_FLAG "--with-precision=single")
else()
    set(PETSC_PRECISION_FLAG "--with-precision=double")
endif()


find_program(SLURM_SRUN_COMMAND srun DOC "Path to the SLURM srun executable")

if (SLURM_SRUN_COMMAND)
    set(PETSC_MPIEXEC ${SLURM_SRUN_COMMAND})
    message(STATUS "Slurm srun detected: ${SLURM_SRUN_COMMAND}")
elseif(${MPIEXEC_EXECUTABLE})
    set(PETSC_MPIEXEC ${MPIEXEC_EXECUTABLE}) # used in cmake > 3.10
    message(STATUS "MPIEXEC_EXECUTABLE detected: ${MPIEXEC_EXECUTABLE}")
elseif(${MPIEXEC})
    set(PETSC_MPIEXEC ${MPIEXEC}) # used in cmake > 3.10
    message(STATUS "MPIEXEC detected: ${MPIEXEC}")
else()
    set(PETSC_MPIEXEC "mpiexec")
endif()


message(STATUS "Slurm executable: ${SLURM_EXECUTABLE}")


if(MADG_USE_SERIAL OR MADG_USE_OPENMP)
ExternalProject_Add(
        petsc_external
        GIT_REPOSITORY https://github.com/petsc/petsc.git
        GIT_TAG 6e6ecd73b34105bf31fbcb05a63437061cbea93b

        BUILD_IN_SOURCE 1
        SOURCE_DIR=${CMAKE_BINARY_DIR}/external/petsc/

        CONFIGURE_COMMAND
        ${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external/configure
        PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external
        PETSC_ARCH=${PETSC_ARCH_FLAG}
        --with-cc=${CMAKE_C_COMPILER} --with-cxx=${CMAKE_CXX_COMPILER} --with-fc=0 --with-pic=1 --with-cxx-dialect=C++17 MAKEFLAGS=$MAKEFLAGS COPTFLAGS=${PETSC_OPT_FLAGS} CXXOPTFLAGS=${PETSC_OPT_FLAGS} --with-mpi=0 --with-debugging=${PETSC_DEBUGGING} --download-hwloc=1 --download-f2cblaslapack=1 ${PETSC_PRECISION_FLAG}

        BUILD_COMMAND
        make -j PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external PETSC_ARCH=${PETSC_ARCH_FLAG}

        INSTALL_COMMAND ""
)

elseif(MADG_USE_CUDA)
    ExternalProject_Add(
            petsc_external
            GIT_REPOSITORY https://github.com/petsc/petsc.git
            GIT_TAG 6e6ecd73b34105bf31fbcb05a63437061cbea93b

            BUILD_IN_SOURCE 1
            SOURCE_DIR=${CMAKE_BINARY_DIR}/external/petsc/

            CONFIGURE_COMMAND
            ${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external/configure
            PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external
            PETSC_ARCH=${PETSC_ARCH_FLAG}
            --with-cc=${CMAKE_C_COMPILER} --with-cxx=${CMAKE_CXX_COMPILER} --with-fc=0 --with-pic=1 --with-cxx-dialect=C++17 MAKEFLAGS=$MAKEFLAGS COPTFLAGS=${PETSC_OPT_FLAGS} CXXOPTFLAGS=${PETSC_OPT_FLAGS} --with-mpi=0 --with-debugging=${PETSC_DEBUGGING} --download-hwloc=1 --download-f2cblaslapack=1 --with-cudac=${CMAKE_CUDA_COMPILER} --with-cuda --with-cuda-arch=35;52;80 --download-kokkos --download-kokkos-kernels --with-kokkos-kernels-tpl=0 --with-cuda-dir=${CUDA_TOOLKIT_ROOT_DIR}  ${PETSC_PRECISION_FLAG}

            BUILD_COMMAND
            make -j PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external PETSC_ARCH=${PETSC_ARCH_FLAG}

            INSTALL_COMMAND ""
    )
elseif(MADG_USE_HIP)
    ExternalProject_Add(
            petsc_external
            GIT_REPOSITORY https://gitlab.com/petsc/petsc.git
            GIT_TAG 6e6ecd73b34105bf31fbcb05a63437061cbea93b

            BUILD_IN_SOURCE 1
            SOURCE_DIR=${CMAKE_BINARY_DIR}/external/petsc/

            CONFIGURE_COMMAND
            ${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external/configure
            PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external
            PETSC_ARCH=${PETSC_ARCH_FLAG}
            --with-cxx=${CMAKE_CXX_COMPILER} --with-fc=0 --with-pic=1 --with-cxx-dialect=C++17 MAKEFLAGS=$MAKEFLAGS COPTFLAGS=${PETSC_OPT_FLAGS} CXXOPTFLAGS=${PETSC_OPT_FLAGS} --with-mpi=0 --with-debugging=${PETSC_DEBUGGING} --download-hwloc=1 --download-f2cblaslapack=1 --with-hip=1 --with-hip-arch=gfx906 --with-hip-dir=/opt/rocm-5.4.3 --download-kokkos --download-kokkos-kernels --with-kokkos-kernels-tpl=0 ${PETSC_PRECISION_FLAG}

            BUILD_COMMAND
            make -j PETSC_DIR=${CMAKE_BINARY_DIR}/external/petsc_external-prefix/src/petsc_external PETSC_ARCH=${PETSC_ARCH_FLAG}

            INSTALL_COMMAND ""
    )
endif()

ExternalProject_Get_Property(petsc_external source_dir)

set(PETSC_INCLUDES
        ${source_dir}/include
        ${source_dir}/${PETSC_ARCH_FLAG}/include
        ${CMAKE_BINARY_DIR}/external/petsc/${PETSC_ARCH_FLAG}/include
        CACHE PATH "")

set(PETSC_LIBRARIES
        ${source_dir}/${PETSC_ARCH_FLAG}/lib/libpetsc.so
        CACHE FILEPATH "")

set(METIS_INCLUDE_DIR ${source_dir}/${PETSC_ARCH_FLAG}/include CACHE PATH "")
set(METIS_LIBRARY ${source_dir}/${PETSC_ARCH_FLAG}/lib/libmetis.so CACHE FILEPATH "")
