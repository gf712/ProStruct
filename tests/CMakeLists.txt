add_subdirectory("${PROJECT_SOURCE_DIR}/third_party/googletest" "extern/googletest")

configure_file(${PROJECT_SOURCE_DIR}/tests/test.pdb ${CMAKE_BINARY_DIR}/tests/test.pdb COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/tests/test.pdb ${CMAKE_BINARY_DIR}/test.pdb COPYONLY)

macro(package_add_test TESTNAME)
    add_executable(${TESTNAME} gtest_suite.cpp ${ARGN})
    target_link_libraries(${TESTNAME} gtest prostruct)
    target_include_directories(${TESTNAME} PRIVATE ${PROJECT_SOURCE_DIR}/src)
    foreach(TEST_FILE ${ARGN})
        get_filename_component(TEST_FILENAME ${TEST_FILE} NAME_WE)
        add_test(NAME ${TEST_FILENAME}
                COMMAND ${CMAKE_BINARY_DIR}/tests/${TESTNAME} --gtest_filter=${TEST_FILENAME}*)
    endforeach()
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

file(GLOB_RECURSE TEST_FILES ${PROJECT_SOURCE_DIR}/tests/*/*.cpp)

package_add_test(prostruct-unit-tests ${TEST_FILES})
