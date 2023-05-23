add_executable(verification_testing
        verification_testing_main.cpp
        poisson/poisson_multigrid_test.cpp
        linear_solver/linear_solver_test.cpp
        )

if(MADG_USE_CUDA)
    message(STATUS "Building poisson/poisson_multigrid_test.cpp with CUDA support.")
    set_source_files_properties(
            poisson/poisson_multigrid_test.cpp
            PROPERTIES LANGUAGE CUDA)
endif()


add_dependencies(verification_testing ${external_libraries})

target_link_libraries(verification_testing
        ${maDG_library_links}
        Boost::program_options
        ${GTEST_LIBRARY_PATH}
        dl
        RAJA
        ${OpenMP_CXX_LIBRARIES}
        )

set_target_properties(verification_testing PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        )

add_test(
        NAME
        unit_testing
        COMMAND
        unit_testing
)

