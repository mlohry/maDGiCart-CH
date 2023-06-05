#########################################################
message( STATUS "Building external Eigen project.")
#########################################################

ExternalProject_Add( eigen_external
#        URL "https://gitlab.com/libeigen/eigen/-/archive/3.3.4/eigen-3.3.4.tar.gz"

        # For CUDA Eigen needed something on newer branches to compile correctly
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 35a367d557078462a0793c88c44dcad64fc63698

        # header only, just extract.
#        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/eigen/
#        UPDATE_COMMAND ""
#        BUILD_COMMAND ${CMAKE_COMMAND} -E copy_directory
#        ${CMAKE_BINARY_DIR}/external/eigen/Eigen
#        ${EXTERNAL_INSTALL_DIR}/include/Eigen
#
#        CONFIGURE_COMMAND ""
#        INSTALL_COMMAND ""
#
#
        CONFIGURE_COMMAND ""
        INSTALL_COMMAND ""
        BUILD_COMMAND ""
        )

#set( EIGEN3_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/eigen PARENT_SCOPE )
set( EIGEN3_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/eigen_external-prefix/src/eigen_external/ PARENT_SCOPE )
