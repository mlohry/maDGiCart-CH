# Get git info and add to generated source
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/git_info.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp" @ONLY)


set(logger_src
        logger.hpp
        logger.cpp
        git_info.hpp
        ${CMAKE_CURRENT_BINARY_DIR}/git_info.cpp
        profiler.cpp profiler.hpp)

add_library(logger ${LIBTYPE} ${logger_src})

add_dependencies(logger ${external_libraries})

target_link_libraries(logger
        Boost::regex
        Boost::log_setup
        Boost::log
        Boost::filesystem
        Boost::thread
        parallel
        )
