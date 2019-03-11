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

This runs very slowly (~500 times slower than the `calculate_phi1` method of `prostruct.PDB`), because numpy is more powerful when it can vectorise calculations, and in this case we are breaking down the problem into N<sub>residues</sub>. The `np.cross` calls are particularly costly! 

Replacing the numpy cross and dot product with a pure python function for 3D vectors will provide speed gains:

```python
import prostruct
import numpy as np
import math

def dot3(a, b):
    return sum([(a_i*b_i) for a_i, b_i in zip(a, b)])

def cross3(a, b):
    return [a[1]*b[2] - a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def dihedral(a1, a2, a3, a4):
    d1 = a1-a2
    d2 = a2-a3
    d3 = a3-a4
    b1 = (d1) / np.linalg.norm(d1)
    b2 = (d2) / np.linalg.norm(d2)
    b3 = (d3) / np.linalg.norm(d3)
    n1 = cross3(b1, b2)
    n2 = cross3(b2, b3)
    return math.degrees(math.atan2(dot3(cross3(n1, b2), n2), dot3(n1, n2)))

class MyStruct(prostruct.CustomPDB): 
    def custom_kernel(self, a, b):
        a = a.get_backbone_atoms()
        b = b.get_backbone_atoms()
        return dihedral(a[2,:], b[0,:], b[1,:], b[2,:])

pdb = MyStruct("mypdb.pdb")
pdb.run_custom_kernel()
```