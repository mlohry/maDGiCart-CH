set(parallel_src
        machine.cpp
        )

if (MADG_USE_CUDA)
    message(STATUS "Building libparallel with CUDA support.")
    set_source_files_properties(
            ${parallel_src}
            PROPERTIES LANGUAGE CUDA)
else()
    message(STATUS "Building libparallel without CUDA support.")
endif()

add_library(parallel ${LIBTYPE} ${parallel_src})

target_link_libraries(parallel
        discretization
        ${OpenMP_CXX_LIBRARIES}
        )

add_dependencies(parallel ${external_libraries})

