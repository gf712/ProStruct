cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 11)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Boost 1.59 REQUIRED COMPONENTS unit_test_framework)

add_subdirectory(prostruct)