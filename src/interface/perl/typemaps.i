%module prostruct
%{
#ifdef SWIGPERL
// name clashes with Perl.h
#ifdef THR
#undef THR
#endif
#ifdef cx_type
#undef cx_type
#endif
#ifdef PP
#undef PP
#endif
#endif
#include "prostruct/prostruct.h"
%}

%define ARMA_COL_OUT(TYPE, PERL_TYPE, PERL_SETTER)
%typemap(out) arma::Col<TYPE>
{
    AV* av = newAV();
    for (int i = 0; i < $1.size(); ++i) 
    {
        SV* perlval = PERL_TYPE(0);
        PERL_SETTER(perlval, $1(i));
        av_push(av, perlval);
    }
    $result = newRV_noinc((SV*) av );
    sv_2mortal( $result );
    argvi++;
}
%enddef

%define ARMA_MAT_OUT(TYPE, PERL_TYPE, PERL_SETTER)
%typemap(out) arma::Mat<TYPE>
{
	AV* av_outer = newAV();    
    for (int i = 0; i < $1.n_cols; ++i) 
    {
    	SV* first = PERL_TYPE(0);
    	AV* av_inner = newAV();

    	for (int j = 0; j < $1.n_rows; ++j) 
    	{
	        SV* perlval = PERL_TYPE(0);
			if (j==0)
    			first = perlval;
	        PERL_SETTER(perlval, $1(i));
	        av_push(av_inner, perlval);
    	}
    	av_push(av_outer, first);
	}
    $result = newRV_noinc((SV*) av_outer );
    sv_2mortal( $result );
    argvi++;
}
%enddef

ARMA_COL_OUT(float, newSVnv, sv_setnv)
ARMA_COL_OUT(double, newSVnv, sv_setnv)

ARMA_MAT_OUT(float, newSVnv, sv_setnv)
ARMA_MAT_OUT(double, newSVnv, sv_setnv)



%include "prostruct.i"