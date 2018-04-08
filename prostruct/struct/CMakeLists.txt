project(struct)

set(TARGET_NAME STRUCT_LIB)

set(SOURCE_FILES
        src/atom.cpp
        src/bond.cpp
        src/residue.cpp
        src/chain.cpp
        )


add_library(${TARGET_NAME} STATIC
        ${SOURCE_FILES}
        )

set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME "struct")

target_include_directories(${TARGET_NAME} SYSTEM PUBLIC
        include/)

target_link_libraries(${TARGET_NAME} PUBLIC
        ${ARMADILLO_LIBRARIES}
        )

add_subdirectory(tests)

add_test(NAME AtomTests COMMAND AtomTest
         WORKING_DIRECTORY ${StructTests_BINARY_DIR})

add_test(NAME BondTests COMMAND BondTest
         WORKING_DIRECTORY ${StructTests_BINARY_DIR})

add_test(NAME ResidueTests COMMAND ResidueTest
        WORKING_DIRECTORY ${StructTests_BINARY_DIR})

add_test(NAME ChainTests COMMAND ChainTest
        WORKING_DIRECTORY ${StructTests_BINARY_DIR})