include( ExternalGTest.cmake )
include( ExternalTabulate.cmake )
include( ExternalRAJA.cmake )

set( external_libraries
        gtest_external
        raja_external
        tabulate_external
        PARENT_SCOPE
        )
