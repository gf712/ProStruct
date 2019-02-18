[![Build Status](https://travis-ci.org/gf712/ProStruct.svg?branch=master)](https://travis-ci.org/gf712/ProStruct)

# ProStruct
ProStruct is a protein structure analysis tool with the end goal of becoming a solution for protein structure analysis
with high performance computing.

To achieve fast speeds and maximise performance given the available linear algebra libraries ProStruct uses
Armadillo. Armadillo can currently use BLAS/LAPACK, OpenBLAS, MKL and NVBLAS. In addition the code is parallelised
where possible using OpenMP (with intra and intercore optimisations).

In addition, ProStruct is available in Python (currently limited) using SWIG. The list of interfaces will
continue to grow.

I also try to write C++ code using the latest standards, so you might need a relatively recent compiler to build
Prostruct.

Contributions are welcome!

## Examples:

C++:
```cpp
#include <prostruct/prostruct.h>

auto pdb = PDB<float>("test.pdb");
auto ks = pdb.compute_kabsch_sander() // arma::Mat<float>
```

Python:
```python
import prostruct

pdb = prostruct.PDB_float("mypdb.pdb")
ks = pdb.compute_kabsch_sander() # numpy array
```

## Build with CMake

Currently ProStruct is only available from source.

```bash
mkdir build
cd build
cmake ..
make
```

With Python interface:
```bash
cmake -DPYTHON_EXECUTABLE=/my/path/to/python -DPYTHON_LIBRARY=/my/path/python/to/lib/libpython3.6m.so -DPYTHON_INCLUDE_DIR=/my/path/to/include/python3.6m/ ..
```