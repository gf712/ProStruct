project(parsers)

set(TARGET_NAME PARSER_LIB)

set(SOURCE_FILES
        src/PDBparser.cpp
        )

add_library(${TARGET_NAME} STATIC
        ${SOURCE_FILES}
        )

set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME "parsers")

target_include_directories(${TARGET_NAME} SYSTEM PUBLIC
        include/)
target_link_libraries(${TARGET_NAME} PUBLIC
        STRUCT_LIB)