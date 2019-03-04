%module prostruct
%{
	#include "prostruct/prostruct.h"
%}

%include "exception.i" 

%exception { 
  try { 
    $action 
  } catch (const std::exception& e) { 
    SWIG_exception(SWIG_RuntimeError, e.what()); 
  } catch (const char* e) { 
    SWIG_exception(SWIG_RuntimeError, e); 
  } catch (const std::string& e) { 
    SWIG_exception(SWIG_RuntimeError, e.c_str()); 
  } 
}

%include "std_string.i"
%include "std_vector.i"

#ifndef SWIGPYTHON
%template() std::vector<std::string>;
#endif

#ifdef SWIGPYTHON
%include <std_shared_ptr.i>
#else 
%include <boost_shared_ptr.i>
#endif

//#ifdef SWIGPYTHON
%define ADD_SHARED_PTR(class_name)
%shared_ptr(class_name)
%enddef
// #else 
// %define ADD_SHARED_PTR(class_name)
// %boost_shared_ptr(class_name)
// %enddef
// #endif

ADD_SHARED_PTR(prostruct::Residue<float>)
ADD_SHARED_PTR(prostruct::Residue<double>)
ADD_SHARED_PTR(prostruct::Chain<float>)
ADD_SHARED_PTR(prostruct::Chain<double>)

//%shared_ptr(prostruct::Residue<float>)
//%shared_ptr(prostruct::Residue<double>)
//%shared_ptr(prostruct::Chain<float>)
//%shared_ptr(prostruct::Chain<double>)

//%shared_ptr(prostruct::PDB)

%include "prostruct/pdb/PDB.h"
%include "prostruct/struct/chain.h"
%include "prostruct/struct/residue.h"

namespace prostruct {
	%template(PDB_float) PDB<float>;
	%template(PDB_double) PDB<double>;

	%template(Chain_float) Chain<float>;
	%template(Chain_double) Chain<double>;

	%template(Residue_float) Residue<float>;
	%template(Residue_double) Residue<double>;
}

%init %{
#ifdef SWIGPYTHON
	import_array();
#endif
%}