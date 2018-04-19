cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -fopenmp-simd -march=native")

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Boost 1.61 REQUIRED COMPONENTS unit_test_framework)
find_package(Armadillo 7.9 REQUIRED)

add_subdirectory(prostruct)