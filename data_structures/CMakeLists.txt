set(data_structures_src
        memory_manager.cpp
        csr_matrix.cpp
        ldu_csr_matrix.cpp
        )

if(MADG_USE_CUDA)
    message(STATUS "Building libdata_structures with CUDA support.")
    set_source_files_properties(
            ${data_structures_src}
            PROPERTIES LANGUAGE CUDA)
else()
    message(STATUS "Building libdata_structures without CUDA support.")
endif()

add_library(data_structures ${LIBTYPE} ${data_structures_src})

add_dependencies(data_structures ${external_libraries})

target_link_libraries(data_structures ${CUDA_LIBRARIES})

if(MADG_USE_HIP)
    target_link_libraries(data_structures roc::rocsparse)

endif()