set(initialization_src
        puppeteer.hpp
        puppeteer.cpp
        )

add_library(initialization ${LIBTYPE} ${initialization_src})

target_link_libraries(initialization ${Boost_LIBRARIES} governing_equations time_stepping discretization fileio ${MATPLOT_LIBRARY})

add_dependencies(initialization ${external_libraries})

