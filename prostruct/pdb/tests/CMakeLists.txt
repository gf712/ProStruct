project(PDBTests)

configure_file(test.pdb ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

add_executable(PDBTest
        PDBTestSuite.cpp)

target_link_libraries(PDBTest
        PDB_LIB
        Boost::unit_test_framework
        )