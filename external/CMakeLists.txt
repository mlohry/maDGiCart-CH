include( ExternalGTest.cmake )
include( ExternalEigen.cmake )
include( ExternalTabulate.cmake )
include( ExternalCelero.cmake )
include( ExternalRAJA.cmake )
include( ExternalMatplotplusplus.cmake )

set( external_libraries
        gtest_external
        eigen_external
        raja_external
        tabulate_external
        celero_external
        matplotplusplus_external
        PARENT_SCOPE
        )
