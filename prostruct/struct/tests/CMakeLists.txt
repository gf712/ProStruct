project(StructTests)

add_executable(AtomTest
        atomTestSuite.cpp)

add_executable(BondTest
        bondTestSuite.cpp)

add_executable(ResidueTest
        residueTestSuite.cpp)

add_executable(ChainTest
        chainTestSuite.cpp)

add_executable(main
        main.cpp)

target_link_libraries(AtomTest
        STRUCT_LIB
        Boost::unit_test_framework
        )

target_link_libraries(BondTest
        STRUCT_LIB
        Boost::unit_test_framework
        )

target_link_libraries(ResidueTest
        STRUCT_LIB
        Boost::unit_test_framework
        )

target_link_libraries(ChainTest
        STRUCT_LIB
        Boost::unit_test_framework
        )

target_link_libraries(main
        STRUCT_LIB
        )