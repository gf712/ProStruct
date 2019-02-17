# - Find the NumPy libraries
# This module finds if NumPy is installed, and sets the following variables
# indicating where it is.
# #
# NUMPY_FOUND - was NumPy found
# NUMPY_VERSION - the version of NumPy found as a string
# NUMPY_VERSION_MAJOR - the major version number of NumPy
# NUMPY_VERSION_MINOR - the minor version number of NumPy
# NUMPY_VERSION_PATCH - the patch version number of NumPy
# NUMPY_VERSION_DECIMAL - e.g. version 1.6.1 is 10601
# NUMPY_INCLUDE_DIRS - path to the NumPy include files

set(NUMPY_FOUND TRUE)

get_filename_component(_python_abs_name "${PYTHON_EXECUTABLE}" ABSOLUTE)

execute_process(COMMAND "${_python_abs_name}" -c
        "import numpy as n; print(n.__version__); print(n.get_include());"
        RESULT_VARIABLE _NUMPY_SEARCH_SUCCESS
        OUTPUT_VARIABLE _NUMPY_VALUES_OUTPUT
        ERROR_VARIABLE _NUMPY_ERROR_VALUE
        OUTPUT_STRIP_TRAILING_WHITESPACE)

if(NOT _NUMPY_SEARCH_SUCCESS MATCHES 0)
    set(NUMPY_FOUND FALSE)
endif()

if (NUMPY_FOUND)
    # Convert the process output into a list
    string(REGEX REPLACE ";" "\\\\;" _NUMPY_VALUES ${_NUMPY_VALUES_OUTPUT})
    string(REGEX REPLACE "\n" ";" _NUMPY_VALUES ${_NUMPY_VALUES})
    list(GET _NUMPY_VALUES 0 NUMPY_VERSION)
    list(GET _NUMPY_VALUES 1 NUMPY_INCLUDE_DIRS)

    # Make sure all directory separators are '/'
    string(REGEX REPLACE "\\\\" "/" NUMPY_INCLUDE_DIRS ${NUMPY_INCLUDE_DIRS})

    # Get the major and minor version numbers
    string(REGEX REPLACE "\\." ";" _NUMPY_VERSION_LIST ${NUMPY_VERSION})
    list(GET _NUMPY_VERSION_LIST 0 NUMPY_VERSION_MAJOR)
    list(GET _NUMPY_VERSION_LIST 1 NUMPY_VERSION_MINOR)
    list(GET _NUMPY_VERSION_LIST 2 NUMPY_VERSION_PATCH)
    string(REGEX MATCH "[0-9]*" NUMPY_VERSION_PATCH ${NUMPY_VERSION_PATCH})
    math(EXPR NUMPY_VERSION_DECIMAL
            "(${NUMPY_VERSION_MAJOR} * 10000) + (${NUMPY_VERSION_MINOR} * 100) + ${NUMPY_VERSION_PATCH}")
endif()

find_package_handle_standard_args(Numpy
        FOUND_VAR NUMPY_FOUND
        REQUIRED_VARS NUMPY_INCLUDE_DIRS
        VERSION_VAR NUMPY_VERSION)

set(NUMPY_FOUND TRUE)