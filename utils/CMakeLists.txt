set(utils_src
        asymptotic_limit.cpp
        )


add_library(utils ${LIBTYPE} ${utils_src})

add_dependencies(utils ${external_libraries})
