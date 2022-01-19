#########################################################
message( STATUS "Building external least-squares-cpp project.")
#########################################################


ExternalProject_Add(
        lsq_external
        URL "https://github.com/Rookfighter/least-squares-cpp/archive/refs/heads/master.zip"
        INSTALL_COMMAND ""

        # header only, just extract.
        SOURCE_DIR ${CMAKE_BINARY_DIR}/external/lsq/
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

set(LSQ_INCLUDE_DIR ${CMAKE_BINARY_DIR}/external/lsq/include/ CACHE PATH "")
