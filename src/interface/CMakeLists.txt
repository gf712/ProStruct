cmake_minimum_required(VERSION 3.13)

find_package(SWIG 3.0 REQUIRED)

set(PYTHON_TARGET_NAME prostruct_python_interface)
set(R_TARGET_NAME prostruct_r_interface)
set(PERL_TARGET_NAME prostruct_perl_interface)


cmake_policy(SET CMP0078 NEW)

if(${CMAKE_VERSION} GREATER "3.14.0")
    cmake_policy(SET CMP0086 OLD)
endif()

include(${SWIG_USE_FILE})
include_directories(swig)
include_directories(${PROJECT_SOURCE_DIR}/src)

#set(CMAKE_SWIG_OUTDIR ${CMAKE_BINARY_DIR}/src/interface/python)
if (ENABLE_PYTHON)
    find_package(PythonLibs 3.6 REQUIRED)
    find_package(Numpy 1.7 REQUIRED)
	file(GLOB_RECURSE SWIG_FILES ${PROJECT_SOURCE_DIR}/src/interface/python/*.i)
	set_property(SOURCE ${SWIG_FILES} PROPERTY CPLUSPLUS ON)

	set(CMAKE_SWIG_FLAGS "-builtin")
	swig_add_library(${PYTHON_TARGET_NAME}
        	LANGUAGE python SOURCES ${SWIG_FILES}
        	OUTFILE_DIR "${CMAKE_BINARY_DIR}/src/interface/python"
        	OUTPUT_DIR "${CMAKE_BINARY_DIR}/src/interface/python")
	target_include_directories(${PYTHON_TARGET_NAME} 
			PUBLIC ${PYTHON_INCLUDE_DIRS} 
			PUBLIC ${PROJECT_SOURCE_DIR}/src)
	if (APPLE)
    	target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC "-undefined dynamic_lookup")
	endif()
	target_link_libraries(${PYTHON_TARGET_NAME} PUBLIC prostruct ${ARMADILLO_LIBRARIES})
	target_include_directories(${PYTHON_TARGET_NAME} PUBLIC ${NUMPY_INCLUDE_DIRS})
	set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES SWIG_DEPENDS prostruct)
	set_target_properties(${PYTHON_TARGET_NAME} PROPERTIES 
						  OUTPUT_NAME prostruct
						  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src/interface/python")
endif()
if (ENABLE_R)
	find_package(R REQUIRED)
	find_package(Boost REQUIRED)
	file(GLOB_RECURSE SWIG_FILES ${PROJECT_SOURCE_DIR}/src/interface/r/*.i)
	set_property(SOURCE ${SWIG_FILES} PROPERTY CPLUSPLUS ON)
	set(CMAKE_SWIG_FLAGS "")
	swig_add_library(${R_TARGET_NAME}
        	LANGUAGE r SOURCES ${SWIG_FILES}
        	OUTFILE_DIR "${CMAKE_BINARY_DIR}/src/interface/r"
        	OUTPUT_DIR "${CMAKE_BINARY_DIR}/src/interface/r")	
	target_include_directories(${R_TARGET_NAME} 
			PUBLIC ${R_INCLUDE_DIR}
			PUBLIC ${PROJECT_SOURCE_DIR}/src
			PUBLIC ${Boost_INCLUDE_DIRS})
	if (APPLE)
    	target_link_libraries(${R_TARGET_NAME} PUBLIC "-undefined dynamic_lookup")
	endif()
	target_link_libraries(${R_TARGET_NAME} PUBLIC prostruct ${ARMADILLO_LIBRARIES} ${R_LIBRARIES})
	set_target_properties(${R_TARGET_NAME} PROPERTIES SWIG_DEPENDS prostruct)
	set_target_properties(${R_TARGET_NAME} PROPERTIES
						  OUTPUT_NAME prostruct
						  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src/interface/r")
endif()

if (ENABLE_PERL)
	find_package(PerlLibs REQUIRED)
	file(GLOB_RECURSE SWIG_FILES ${PROJECT_SOURCE_DIR}/src/interface/perl/*.i)
	set_property(SOURCE ${SWIG_FILES} PROPERTY CPLUSPLUS ON)
	set(CMAKE_SWIG_FLAGS "")
	swig_add_library(${PERL_TARGET_NAME}
        	LANGUAGE Perl SOURCES ${SWIG_FILES}
        	OUTFILE_DIR "${CMAKE_BINARY_DIR}/src/interface/perl"
        	OUTPUT_DIR "${CMAKE_BINARY_DIR}/src/interface/perl")	
	target_include_directories(${PERL_TARGET_NAME} 
			PUBLIC ${PERL_INCLUDE_PATH}
			PUBLIC ${PROJECT_SOURCE_DIR}/src)
	if (APPLE)
    	target_link_libraries(${PERL_TARGET_NAME} PUBLIC "-undefined dynamic_lookup")
	endif()
	target_link_libraries(${PERL_TARGET_NAME} PUBLIC prostruct ${ARMADILLO_LIBRARIES} ${PERL_LIBRARY} "-flat_namespace")
	set_target_properties(${PERL_TARGET_NAME} PROPERTIES SWIG_DEPENDS prostruct)
	set_target_properties(${PERL_TARGET_NAME} PROPERTIES
						  OUTPUT_NAME prostruct
						  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/src/interface/perl")
endif()