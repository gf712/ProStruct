set(TARGET_NAME prostruct)

add_subdirectory("${PROJECT_SOURCE_DIR}/third_party/fmt" "extern/fmt")

file(GLOB_RECURSE SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/prostruct/*.cpp)

add_library(${TARGET_NAME} SHARED
        ${SOURCE_FILES}
        )

set(PROSTRUCT_INCLUDE ${PROJECT_SOURCE_DIR}/src)

target_include_directories(${TARGET_NAME}
        PUBLIC ${PROSTRUCT_INCLUDE}
        PUBLIC ${CMAKE_BINARY_DIR}/src)

target_link_libraries(${TARGET_NAME}
        PUBLIC ${ARMADILLO_LIBRARIES}
        PUBLIC fmt::fmt)