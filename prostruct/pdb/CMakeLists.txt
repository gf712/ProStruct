project(pdb)

set(TARGET_NAME PDB_LIB)

set(SOURCE_FILES
        src/PDB.cpp
        src/dssp.cpp
        )

add_library(${TARGET_NAME} STATIC
        ${SOURCE_FILES}
        )

set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME "pdb")

target_include_directories(${TARGET_NAME}
        INTERFACE include/)

target_link_libraries(${TARGET_NAME}
        PARSER_LIB
        STRUCT_LIB
        )

add_subdirectory(tests)

add_test(NAME PDBTests COMMAND PDBTest
        WORKING_DIRECTORY ${PDBTests_BINARY_DIR})