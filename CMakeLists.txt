cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-O3 -march=native") # -fopenmp -fopenmp-simd ")

option(ENABLE_TESTS "Enable tests" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Armadillo 7.9 REQUIRED)

add_subdirectory(src/prostruct)

IF(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
ENDIF()

configure_file(${PROJECT_SOURCE_DIR}/tests/test.pdb ${CMAKE_BINARY_DIR}/test.pdb COPYONLY)
add_executable(load_test tests/main.cpp)
target_include_directories(load_test PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(load_test prostruct)