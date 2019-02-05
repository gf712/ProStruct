cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fopenmp-simd -march=native")

option(ENABLE_TESTS "Enable tests" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Armadillo 7.9 REQUIRED)

add_subdirectory(prostruct)

IF(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
ENDIF()