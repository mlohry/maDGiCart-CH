include( ExternalGTest.cmake )
include( ExternalTabulate.cmake )
include( ExternalRAJA.cmake )
include( ExternalMatplotplusplus.cmake )

set( external_libraries
        gtest_external
        raja_external
        tabulate_external
        matplotplusplus_external
        PARENT_SCOPE
        )
