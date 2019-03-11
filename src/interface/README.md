# Extending the Python functionality

At its core most of the calculations performed in ProStruct go via the [kernel engine](../prostruct/README.md). To further extend the functionality ProStruct offers and extendable class called CustomPDB which is available in the Python interface. It is a SWIG director which allows overriding base class methods in the interface language. In this case the user can write a new `custom_kernel` method that returns a scalar which is then run using `run_custom_kernel`.

Here is an example of how to calculate the phi angles with SWIG directors:

```python
import prostruct
import numpy as np

def dihedral(a1, a2, a3, a4):
	b1 = (a1-a2) / np.linalg.norm(a1-a2)
	b2 = (a2-a3) / np.linalg.norm(a2-a3)
	b3 = (a3-a4) / np.linalg.norm(a3-a4)
	n1 = np.cross(b1, b2)
	n2 = np.cross(b2, b3)
	result = np.rad2deg(np.arctan2(np.dot(np.cross(n1, b2), n2), np.dot(n1, n2)))
	return result

class MyStruct(prostruct.CustomPDB): 
   def custom_kernel(self, a, b):
   		a = a.get_backbone_atoms()
   		b = b.get_backbone_atoms()
   		return float(dihedral(a[2,:], b[0,:], b[1,:], b[2,:]))

pdb = MyStruct("mypdb.pdb")
pdb.run_custom_kernel()
```

Currently this runs very slowly (~500 times slower than the `calculate_phi1` method of PDB), because there are lots of conversions between the C++ and Python objects, i.e. armadillo matrices to python's numpy arrays.