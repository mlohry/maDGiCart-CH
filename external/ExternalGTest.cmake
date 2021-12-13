#########################################################
message(STATUS "Building external GTest project.")
#########################################################

ExternalProject_Add(
        gtest_external
        URL "https://github.com/google/googletest/archive/release-1.8.0.zip"
        URL_MD5 adfafc8512ab65fd3cf7955ef0100ff5
        INSTALL_COMMAND ""
)

# Specify include dir
ExternalProject_Get_Property(gtest_external source_dir)

set(GTEST_INCLUDE_DIR
        ${source_dir}/googletest/include
        CACHE PATH "")

# Library
ExternalProject_Get_Property(gtest_external binary_dir)

set(GTEST_LIBRARY_PATH
        ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest.a
        CACHE FILEPATH "")

set(GTEST_LIBRARY_MAIN_PATH
        ${binary_dir}/googlemock/gtest/${CMAKE_FIND_LIBRARY_PREFIXES}gtest_main.a
        CACHE FILEPATH "")

set(GTEST_LIBRARY gtest)

add_library(${GTEST_LIBRARY} UNKNOWN IMPORTED)

set_property(TARGET ${GTEST_LIBRARY} PROPERTY IMPORTED_LOCATION
        ${GTEST_LIBRARY_PATH})

add_dependencies(${GTEST_LIBRARY} gtest_external)
