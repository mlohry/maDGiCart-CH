#########################################################
message( STATUS "Building external Matplot++ project.")
#########################################################


ExternalProject_Add( matplotplusplus_external

        URL "https://github.com/alandefreitas/matplotplusplus/archive/refs/tags/v1.1.0.zip"
        URL_MD5 10c6161f20af9fa556ecd120a5db708d

        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/matplotplusplus
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/matplotplusplus-build

        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_DIR} -DBUILD_SHARED_LIBS=On -DBUILD_TESTS=Off -DBUILD_TESTING=Off -DBUILD_EXAMPLES=Off

        )


set(MATPLOT_LIBRARY
        ${EXTERNAL_INSTALL_DIR}/lib/libmatplot.so
        CACHE PATH ""
        )
