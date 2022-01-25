#########################################################
message( STATUS "Building external TinyXML2 project.")
#########################################################


ExternalProject_Add(
        tinyxml2_external
        URL "https://github.com/leethomason/tinyxml2/archive/refs/tags/9.0.0.zip"
        URL_MD5 2a3b1b8acdc1a0bd15e4010d91c505f8

        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/tinyxml2/
        BINARY_DIR ${CMAKE_BINARY_DIR}/external/tinyxml2-build/
        INSTALL_DIR ${EXTERNAL_INSTALL_DIR}
        CMAKE_ARGS -D CMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_DIR} -D tinyxml2_BUILD_TESTING=Off -DBUILD_SHARED_LIBS=On
)

ExternalProject_Get_Property(tinyxml2_external source_dir)

set(TINYXML2_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/tinyxml2/ CACHE PATH "")
set(TINYXML2_LIBRARIES ${CMAKE_BINARY_DIR}/external/tinyxml2-build/libtinyxml2.so CACHE PATH "")
