set(governing_equations_src
        cahn_hilliard/cahn_hilliard_2d_fd.cpp
        cahn_hilliard/cahn_hilliard_initial_conditions.cpp
        residual_calculator.cpp
        cahn_hilliard/cahn_hilliard_solution_monitor.cpp)


if(MADG_USE_CUDA)
        message(STATUS "Building libgoverning_equations with CUDA support.")
        set_source_files_properties(
                ${governing_equations_src}
         PROPERTIES
                LANGUAGE CUDA
                CUDA_STANDARD 11
        )
         
else()
        message(STATUS "Building libgoverning_equations without CUDA support.")
endif()



add_library(governing_equations ${LIBTYPE} ${governing_equations_src})


add_dependencies(governing_equations ${external_libraries})



if(MADG_USE_CUDA)
        target_link_libraries(governing_equations
                discretization
                parallel
                RAJA
                ${CUDA_LIBRARIES}
                )
else()
        target_link_libraries(governing_equations
                discretization
                parallel
                )
endif()
