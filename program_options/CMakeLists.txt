set(program_options_src
        program_options.hpp
        program_options_parser.cpp
        program_options_parser.hpp
        )

add_library(program_options
        ${LIBTYPE} ${program_options_src}
        )

target_link_libraries(program_options Boost::program_options)

add_dependencies(program_options ${external_libraries})
