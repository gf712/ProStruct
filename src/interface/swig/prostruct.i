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
#elif SWIGR
%include <boost_shared_ptr.i>
#endif

#ifndef SWIGPERL
%shared_ptr(prostruct::Residue<float>)
%shared_ptr(prostruct::Residue<double>)
%shared_ptr(prostruct::Chain<float>)
%shared_ptr(prostruct::Chain<double>)
#endif
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