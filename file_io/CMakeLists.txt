set(fileio_src
        vtk_solution_writer.cpp
        )

add_library(fileio ${LIBTYPE} ${fileio_src})

target_link_libraries(fileio logger)

add_dependencies(fileio ${external_libraries})
