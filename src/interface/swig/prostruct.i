/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */
 
#ifdef SWIGPYTHON
%module(directors="1") prostruct
#else
%module prostruct
#endif

#ifdef SWIGPYTHON
%include <std_shared_ptr.i>
#elif SWIGR
%include <boost_shared_ptr.i>
#endif

%include "exception.i" 
%include "std_string.i"
%include "std_vector.i"

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

#ifndef SWIGPYTHON
%template() std::vector<std::string>;
#endif

%{
  #include <prostruct/prostruct.h>
%}

%include "prostruct/struct/utils.h"
%include "prostruct/struct/chain.h"
%include "prostruct/struct/residue.h"
%include "prostruct/pdb/struct_base.h"
%include "prostruct/struct/chain.h"
%include "prostruct/pdb/PDB.h"

#ifndef SWIGPERL
%shared_ptr(prostruct::StructBase<float>)
%shared_ptr(prostruct::StructBase<double>)
%shared_ptr(prostruct::PDB<float>)
%shared_ptr(prostruct::PDB<double>)
%shared_ptr(prostruct::Chain<float>)
%shared_ptr(prostruct::Chain<double>)
%shared_ptr(prostruct::Residue<float>)
%shared_ptr(prostruct::Residue<double>)
#endif

// TODO: add support for vectors
// %template(VectorOfResidues_float)  std::vector<std::shared_ptr<prostruct::Residue<float>>>;
// %template(VectorOfResidues_double) std::vector<std::shared_ptr<prostruct::Residue<double>>>;
// %template(VectorOfAtoms_float)     std::vector<std::shared_ptr<prostruct::Atom<float>>>;
// %template(VectorOfAtoms_double)    std::vector<std::shared_ptr<prostruct::Atom<double>>>;
// %template(VectorOfChains_float)    std::vector<std::shared_ptr<prostruct::Chain<float>>>;
// %template(VectorOfChains_double)   std::vector<std::shared_ptr<prostruct::Chain<double>>>;

%template(StructBase_float)  prostruct::StructBase<float>;
%template(StructBase_double) prostruct::StructBase<double>;  

%template(PDB_float)  prostruct::PDB<float>;
%template(PDB_double) prostruct::PDB<double>;

%template(Chain_float)  prostruct::Chain<float>;
%template(Chain_double) prostruct::Chain<double>;

%template(Residue_float)  prostruct::Residue<float>;
%template(Residue_double) prostruct::Residue<double>;

#ifdef SWIGPYTHON
%shared_ptr(prostruct::CustomPDB)
%director prostruct::CustomPDB;
%feature("director:except") {
    if ($error != NULL) {
        throw Swig::DirectorMethodException();
    }
}
%include "prostruct/pdb/custom_pdb.h"
#endif

%init %{
#ifdef SWIGPYTHON
	import_array();
#endif
%}