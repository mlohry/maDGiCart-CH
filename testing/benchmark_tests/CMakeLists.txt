add_executable(benchmark_testing
        benchmark_main.cpp
        cahn_hilliard_3d_benchmark.cpp
        cahn_hilliard_2d_benchmark.cpp
        )


add_dependencies(benchmark_testing ${external_libraries})

target_include_directories(benchmark_testing PUBLIC ${CELERO_INCLUDE_DIR})

target_link_libraries(benchmark_testing
        ${maDG_library_links}
        Boost::program_options
        ${GTEST_LIBRARY_PATH}
        dl
        RAJA
        ${OpenMP_CXX_LIBRARIES}
        ${CELERO_LIBRARY}
        )

set_target_properties(benchmark_testing PROPERTIES
       RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
       )