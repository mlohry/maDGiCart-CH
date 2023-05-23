#########################################################
message(STATUS "Building external Celero project.")
#########################################################

ExternalProject_Add(
        celero_external
        URL https://github.com/DigitalInBlue/Celero/archive/refs/tags/v2.8.5.zip
        URL_MD5 50b779625cf9f4ce2206b06859b3182b
        CMAKE_ARGS -DCMAKE_CXX_STANDARD=14
        INSTALL_COMMAND ""
)

set(CELERO_INCLUDE_DIR
        ${CMAKE_BINARY_DIR}/external/celero_external-prefix/src/celero_external/include/
        CACHE PATH "")

message(STATUS "Celero external include dir: ${CELERO_INCLUDE_DIR}")


ExternalProject_Get_Property(celero_external binary_dir)

set(CELERO_LIBRARY
        ${binary_dir}/libcelero.so
        CACHE PATH ""
        )

message(STATUS "Celero library: ${CELERO_LIBRARY}")
