%module prostruct
%{
    #include "prostruct/prostruct.h"
%}

%include "std_string.i"
%include "std_vector.i"

#ifndef SWIGPYTHON
%template() std::vector<std::string>;
#endif

%include <std_shared_ptr.i>

%shared_ptr(prostruct::Chain<float>)
%shared_ptr(prostruct::Chain<double>)
//%shared_ptr(prostruct::PDB)

%include "prostruct/pdb/PDB.h"
%include "prostruct/struct/chain.h"

%template(PDB_float) prostruct::PDB<float>;
%template(PDB_double) prostruct::PDB<double>;

%template(Chain_float) prostruct::Chain<float>;
%template(Chain_double) prostruct::Chain<double>;

%init %{
#ifdef SWIGPYTHON
	import_array();
#endif
%}