set(fileio_src
        vtk_solution_writer.cpp
        vtk_solution_reader.cpp
        )

add_library(fileio ${LIBTYPE} ${fileio_src})

target_link_libraries(fileio logger "${TINYXML2_LIBRARIES}")

add_dependencies(fileio ${external_libraries})
