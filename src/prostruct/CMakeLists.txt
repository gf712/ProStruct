set(TARGET_NAME prostruct)

file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/prostruct/*.cpp)

add_library(${TARGET_NAME} SHARED
        ${SOURCE_FILES}
        )

target_include_directories(${TARGET_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(${TARGET_NAME} PUBLIC
        ${ARMADILLO_LIBRARIES}
        )