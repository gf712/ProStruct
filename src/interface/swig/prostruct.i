%module prostruct
%{
    #include "prostruct/pdb/PDB.h"
%}

%include "std_string.i"
%include "std_vector.i"
%include "prostruct/pdb/PDB.h"

#ifndef SWIGPYTHON
%template() std::vector<std::string>;
#endif

%template(PDB_float) PDB<float>;
%template(PDB_double) PDB<double>;

%init %{
#ifdef SWIGPYTHON
	import_array();
#endif
%}