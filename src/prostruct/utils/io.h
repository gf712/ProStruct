//
// Created by Gil Hoben on 2019-02-18.
//

#ifndef PROSTRUCT_IO_H
#define PROSTRUCT_IO_H

#include "prostruct/config.h"

#ifdef HAVE_CXA_DEMANGLE
#include <cxxabi.h>
#endif

namespace prostruct {
    template <typename T>
    std::string demangled_type() {
        const char* name = typeid(T).name();
#ifdef HAVE_CXA_DEMANGLE
        size_t length;
        int status;
        char *demangled = abi::__cxa_demangle(name, nullptr, &length, &status);
        std::string demangled_string(demangled);
        free(demangled);
#else
        std::string demangled_string(name);
#endif
        return demangled_string;
    }
}

#endif //PROSTRUCT_IO_H
