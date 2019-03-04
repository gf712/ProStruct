%{
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "numpy/ndarrayobject.h"
%}

%feature("python:slot", "tp_repr", functype="reprfunc") prostruct::PDB::to_string;

%typemap(out) std::vector<std::string>
{
	$result = PyList_New($1.size());
	int iter = 0;
	for(auto const& each: *&$1) {
		PyList_SetItem($result, iter, PyUnicode_FromString(each.c_str()));
		iter++;
	}
}

%define ARMA_COL_OUT(TYPE, NUMPY_TYPE)
%typemap(out) arma::Col<TYPE>
{
	npy_intp size[1] = {static_cast<npy_intp>($1.size())};
	$result = PyArray_SimpleNew(1, size, NUMPY_TYPE);
	memcpy(PyArray_DATA((PyArrayObject *) $result), $1.memptr(), sizeof(TYPE) * $1.size());
}
%enddef

%define ARMA_MAT_OUT(TYPE, NUMPY_TYPE)
%typemap(out) arma::Mat<TYPE>
{
	npy_intp size[2] = {static_cast<npy_intp>($1.n_cols), static_cast<npy_intp>($1.n_rows)};
	$result = PyArray_SimpleNew(2, size, NUMPY_TYPE);
	memcpy(PyArray_DATA((PyArrayObject *) $result), $1.memptr(), sizeof(TYPE) * $1.size());
}
%enddef

ARMA_COL_OUT(float, NPY_FLOAT);
ARMA_COL_OUT(double, NPY_DOUBLE);
ARMA_COL_OUT(arma::uword, NPY_ULONGLONG);
ARMA_COL_OUT(unsigned long long, NPY_ULONGLONG);


ARMA_MAT_OUT(float, NPY_FLOAT);
ARMA_MAT_OUT(double, NPY_DOUBLE);
ARMA_MAT_OUT(arma::uword, NPY_ULONGLONG);

%include "prostruct.i"