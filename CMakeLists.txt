cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "-O3 -march=native") # -fopenmp -fopenmp-simd ")

option(ENABLE_TESTS "Enable tests" OFF)
option(ENABLE_INTERFACE "Enable interface" OFF)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Armadillo 7.9 REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_subdirectory(src/prostruct)

IF(ENABLE_INTERFACE)
    find_package(SWIG 3.0 REQUIRED)
    find_package(PythonLibs 3.6 REQUIRED)
    find_package(Numpy 1.7 REQUIRED)
    add_subdirectory(src/interface)
ENDIF()

IF(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(tests)
ENDIF()

configure_file(${PROJECT_SOURCE_DIR}/tests/test.pdb ${CMAKE_BINARY_DIR}/test.pdb COPYONLY)
add_executable(load_test tests/main.cpp)
target_include_directories(load_test PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(load_test prostruct)