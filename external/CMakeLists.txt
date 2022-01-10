include( ExternalGTest.cmake )
include( ExternalTabulate.cmake )
include( ExternalRAJA.cmake )
include( ExternalCelero.cmake )

set( external_libraries
        gtest_external
        raja_external
        tabulate_external
        celero_external
        PARENT_SCOPE
        )
