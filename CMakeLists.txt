cmake_minimum_required(VERSION 3.9)
project(ProStruct)

set(CMAKE_CXX_STANDARD 17)

option(ENABLE_TESTING "Enable tests" OFF)
option(OPTIMIZE_FOR_NATIVE "Build with -march=native" ON)

option(ENABLE_PYTHON "Enable Python interface" OFF)
option(ENABLE_R "Enable R interface" OFF)
option(ENABLE_PERL "Enable Perl interface" OFF)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

find_package(Armadillo 7.9 REQUIRED)

option(USE_OPENMP "Enable OpenMP" ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND (NOT APPLE) AND USE_OPENMP)
    set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
    set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
    set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
    set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
endif()

if ((NOT APPLE) AND USE_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND AND USE_OPENMP)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
elseif(APPLE AND USE_OPENMP)
    MESSAGE(STATUS "OpenMP not supported")
    set(USE_OPENMP OFF)
elseif(NOT USE_OPENMP)
    MESSAGE(STATUS "OpenMP disabled")
endif()

# from https://stackoverflow.com/a/34317588
find_program(CCACHE_EXECUTABLE ccache)
mark_as_advanced(CCACHE_EXECUTABLE)
if(CCACHE_EXECUTABLE)
    foreach(LANG C CXX)
        if(NOT DEFINED CMAKE_${LANG}_COMPILER_LAUNCHER AND NOT CMAKE_${LANG}_COMPILER MATCHES ".*/ccache$")
            message(STATUS "Enabling ccache for ${LANG}")
            set(CMAKE_${LANG}_COMPILER_LAUNCHER ${CCACHE_EXECUTABLE} CACHE STRING "")
        endif()
    endforeach()
endif()

include(CheckCompilerFeatures)

add_subdirectory(src/prostruct)

IF(${ENABLE_PYTHON} OR ${ENABLE_R})
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