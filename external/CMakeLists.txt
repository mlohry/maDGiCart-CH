include( ExternalGTest.cmake )
include( ExternalTabulate.cmake )
include( ExternalRAJA.cmake )
include( ExternalCelero.cmake )
include( ExternalTinyXML2.cmake )

set( external_libraries
        gtest_external
        raja_external
        tabulate_external
        celero_external
        tinyxml2_external
        PARENT_SCOPE
        )
