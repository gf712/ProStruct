include(CheckCXXSourceCompiles)
include(CheckCXXCompilerFlag)

check_cxx_source_compiles(
    "#include <cxxabi.h>
    int main(int argc, char* argv[])
        { char * type; int status;
          char * r = abi::__cxa_demangle(type, 0, 0, &status);
          return 0;
        }"
    HAVE_CXA_DEMANGLE
)

if(HAVE_CXA_DEMANGLE)
    file(WRITE ${CMAKE_BINARY_DIR}/src/prostruct/config.h "#define HAVE_CXA_DEMANGLE 1")
endif()

check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE AND OPTIMIZE_FOR_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()