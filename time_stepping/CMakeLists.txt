set(time_stepping_src
        time_integrator.cpp
        rk3ssp.cpp
        ode23.cpp
        petsc_time_integrator.cpp
        petsc_environment.cpp
        bdf_multigrid.cpp
        ode23_smoother.cpp)

if(MADG_USE_CUDA)
    message(STATUS "Building libtime_stepping with CUDA support.")
    set_source_files_properties(
            ${time_stepping_src}
            PROPERTIES LANGUAGE CUDA)
else()
    message(STATUS "Building libtime_stepping without CUDA support.")
endif()

add_library(time_stepping ${LIBTYPE} ${time_stepping_src})

add_dependencies(time_stepping ${external_libraries})

target_link_libraries(time_stepping parallel logger ${PETSC_LIBRARIES} )

