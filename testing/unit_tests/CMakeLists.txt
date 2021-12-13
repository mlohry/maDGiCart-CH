add_executable(unit_testing
        unit_testing_main.cpp
        spatial_discretization/discretization_2d_cart_test.cpp
        governing_equations/cahn_hilliard_2d_fd_test.cpp
        time_stepping/cahn_hilliard_regression_test.cpp
        )

add_dependencies(unit_testing ${external_libraries})

target_link_libraries(unit_testing
        ${maDG_library_links}
        Boost::program_options
        ${GTEST_LIBRARY_PATH}
        dl
        RAJA
        ${OpenMP_CXX_LIBRARIES}
        )

set_target_properties(unit_testing PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
        )

add_test(
        NAME
        unit_testing
        COMMAND
        unit_testing
)

