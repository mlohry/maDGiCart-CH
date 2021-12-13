#########################################################
message( STATUS "Building external tabulate project.")
#########################################################


ExternalProject_Add(
        tabulate_external
        URL "https://github.com/p-ranav/tabulate/archive/v1.3.zip"
        URL_MD5 743801dc7d47b1c397b279f0eddf257d
        INSTALL_COMMAND ""

        # header only, just extract.
        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/tabulate/
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

set(TABULATE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/tabulate/include/ CACHE PATH "")
