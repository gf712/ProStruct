cmake_minimum_required(VERSION 3.13)

set(TARGET_NAME prostruct_interface)
cmake_policy(SET CMP0078 NEW)

if(${CMAKE_VERSION} GREATER "3.14.0")
    cmake_policy(SET CMP0086 NEW)
endif()


include(${SWIG_USE_FILE})

file(GLOB_RECURSE SWIG_FILES ${PROJECT_SOURCE_DIR}/src/interface/python/*.i)
set_property(SOURCE ${SWIG_FILES} PROPERTY CPLUSPLUS ON)
include_directories(swig)
include_directories(${PROJECT_SOURCE_DIR}/src)

swig_add_library(${TARGET_NAME}
        LANGUAGE python SOURCES ${SWIG_FILES})
target_include_directories(${TARGET_NAME} PUBLIC ${PYTHON_INCLUDE_DIRS} PUBLIC ${PROJECT_SOURCE_DIR}/src)
if (APPLE)
    target_link_libraries(${TARGET_NAME} PUBLIC "-undefined dynamic_lookup")
endif()
target_link_libraries(${TARGET_NAME} PUBLIC prostruct ${ARMADILLO_LIBRARIES})
target_include_directories(${TARGET_NAME} PUBLIC ${NUMPY_INCLUDE_DIRS})
set_target_properties(${TARGET_NAME} PROPERTIES SWIG_DEPENDS prostruct)
set_target_properties(${TARGET_NAME} PROPERTIES OUTPUT_NAME prostruct)