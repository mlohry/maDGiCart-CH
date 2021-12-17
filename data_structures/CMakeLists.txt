set(data_structures_src
        memory_manager.cpp
        )

add_library(data_structures ${LIBTYPE} ${data_structures_src})

add_dependencies(data_structures ${external_libraries})

target_link_libraries(data_structures ${CUDA_LIBRARIES} ${CUDA_CUBLAS_LIBRARIES})
