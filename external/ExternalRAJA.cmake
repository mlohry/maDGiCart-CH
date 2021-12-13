#########################################################
message( STATUS "Building external RAJA project.")
#########################################################


ExternalProject_Add(
        raja_external
        URL https://github.com/LLNL/RAJA/releases/download/v0.14.0/RAJA-v0.14.0.tar.gz
        URL_MD5 1e8682504f0c44301f56c116c14487c2

        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/raja/
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/raja-build/
        INSTALL_DIR ${EXTERNAL_INSTALL_DIR}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_DIR} -DENABLE_CUDA=${MADG_USE_CUDA} -DENABLE_OPENMP=${MADG_USE_OPENMP} -DCMAKE_CUDA_STANDARD=14 -DRAJA_ENABLE_TESTS=Off -DRAJA_ENABLE_EXAMPLES=Off -DRAJA_ENABLE_EXERCISES=Off
)

set(RAJA_LIBRARIES
        ${CMAKE_BINARY_DIR}/external_install/lib/libRAJA.a
        CACHE FILEPATH "")
