/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */
 
%{
	#include <boost/smart_ptr/shared_ptr.hpp>
%}



%define ARMA_COL_OUT(TYPE, R_TYPESXP, R_CAST_TYPE, R_TYPE, R_TYPE_STRING)
%typemap(out) arma::Col<TYPE>
{
   	Rf_protect($result = Rf_allocVector(R_TYPESXP, $1.size()));  
   	printf("TEST");
   	fflush(stdout);
   	for (arma::uword i=0; i< $1.size(); ++i) 
   	{     
   		R_TYPE($result)[i] = static_cast<R_CAST_TYPE>($1(i));
	}
	Rf_unprotect(1);
}

%typemap("rtype") arma::Col<TYPE> R_TYPE_STRING

%typemap("scoerceout") arma::Col<TYPE>
%{ %}

%enddef

%define ARMA_MAT_OUT(TYPE, R_TYPESXP, R_CAST_TYPE, R_TYPE)
%typemap(out) arma::Mat<TYPE>
{
   	Rf_protect($result = Rf_allocMatrix(R_TYPESXP, $1.n_cols, $1.n_rows));    
   	for (arma::uword i=0; i < $1.n_rows; ++i) 
   	{
   		for (arma::uword j=0; j < $1.n_cols; ++j) 
   		{
   			R_TYPE($result)[i*$1.n_cols+j] = static_cast<R_CAST_TYPE>($1(i, j));   
   		}
	}
	Rf_unprotect(1);
}

%typemap("rtype") arma::Mat<TYPE> "matrix"

%typemap("scoerceout") arma::Mat<TYPE>
%{ %}

%enddef

ARMA_COL_OUT(float, REALSXP, float, REAL, "numeric");
ARMA_COL_OUT(double, REALSXP, float, REAL, "numeric");
ARMA_COL_OUT(arma::uword, INTSXP, int, INTEGER, "integer");

ARMA_MAT_OUT(float, REALSXP, float, REAL);
ARMA_MAT_OUT(double, REALSXP, float, REAL);
ARMA_MAT_OUT(arma::uword, INTSXP, int, INTEGER);

%include "prostruct.i"