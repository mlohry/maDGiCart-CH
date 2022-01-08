set(discretization_src
        discretization_2d_cart.cpp
        discretization_3d_cart.cpp)



if(MADG_USE_CUDA)
    message(STATUS "Building libdiscretization with CUDA support.")
    set_source_files_properties(
            ${discretization_src}
            PROPERTIES LANGUAGE CUDA)
else()
    message(STATUS "Building libdiscretization without CUDA support.")
endif()


add_library(discretization ${LIBTYPE} ${discretization_src})

add_dependencies(discretization ${external_libraries})

target_link_libraries(discretization data_structures)

