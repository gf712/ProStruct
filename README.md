[![Build Status](https://travis-ci.org/gf712/ProStruct.svg?branch=master)](https://travis-ci.org/gf712/ProStruct)

# ProStruct
ProStruct is a protein structure analysis tool with the end goal of becoming a solution for protein structure analysis
with high performance computing.

To achieve fast speeds and maximise performance given the available linear algebra libraries ProStruct uses
Armadillo. Armadillo can currently use BLAS/LAPACK, OpenBLAS, MKL and NVBLAS. In addition the code is parallelised
where possible using OpenMP (with intra and intercore optimisations).

In addition, ProStruct is available in Python (currently limited) using SWIG. Hopefully the list of interfaces will
continue to grow.

I also try to write C++ code using the latest standards, so you might need a relatively recent compiler to build
Prostruct.

Contributions are welcome!