
# PDB

This needs to be renamed to something else, i.e. topology.

## How it works:
Currently most of the library is written in (sort of) old school C++ and written like a researcher would.
This is fine, but becomes difficult to expand the library. All the code that supports the PDB class is in `prostruct/struct` and the parser in `struct/parsers`.
Therefore I started working on a friendlier library API, that requires fewer interations with linear algebra code when extending the PDB (to be renamed) class:
- the computation engine:
  Some computations, such as dihedral angle calcuations, only require a residue and its neighbouring residues.
  So this problem can be abstracted to a function that slides along the protein sequence and performs some calculation for each residue(s), i.e. a kernel. 
  This is what the computation engine does: it takes a kernel and slides it along the residue vector.
  The kernel is a C++ lambda that receives Residue object(s) as parameter(s) and returns a scalar value.
  For example, the phi dihedral angle:
  	The computation requires the C atom coordinates of the current residue and the coordinates for the N, CA and C in the following residue. The four atoms form two planes (i.e. C<sub>i</sub>-N<sub>i+1</sub>-CA<sub>i+1</sub> and N<sub>i+1</sub>-CA<sub>i+1</sub>-C<sub>i+1</sub>) and the angle of the intersection is the dihedral angle ([see](https://en.wikipedia.org/wiki/Dihedral_angle)).
  	In C++:
  	```cpp
	auto phi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue, const std::shared_ptr<Residue<T>>& residue_next) {
		if (residue->is_c_terminus())
			return static_cast<T>(0.0);
		// backbone atoms are in this order: N | CA | C | O
		auto atom_coords_this = residue->get_backbone_atoms();
		auto atom_coords_next = residue_next->get_backbone_atoms();
		arma::Col<T2> b1 = arma::normalise(atom_coords_this.col(2) - atom_coords_next.col(0));
		arma::Col<T2> b2 = arma::normalise(atom_coords_next.col(0) - atom_coords_next.col(1));
		arma::Col<T2> b3 = arma::normalise(atom_coords_next.col(1) - atom_coords_next.col(2));
		arma::Col<T2> n1 = arma::cross(b1, b2); // plane 1
		arma::Col<T2> n2 = arma::cross(b2, b3); // plane 2
		// coef is just the factor needed to convert from radians to degrees
		return std::atan2(
				   arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2)) * coef;
	};
  	```
  	The kernel can then be passed to the engine which slides along the sequence:
  	```cpp
  	geometry::residue_kernel_engine(m_residues, 0, phi_kernel);
  	```
  	The compiler determines at compile time the window size required given the number of arguments of the lambda, i.e. 2.
  	The 0 just indicates that we want to start at position 0 in the sequence (N-terminus).
  	The result is a matrix 1 x N<sub>residues</sub>, i.e. a row vector.
  	`geometry::residue_kernel_engine` can take an arbitrary number of kernels. The following would return a 2 x N<sub>residues</sub> where the first row has the phi angles and the second psi angles.
  	```cpp
  	geometry::residue_kernel_engine(m_residues, 0, phi_kernel, psi_kernel);`
  	```
  	The advantage of this approach is that the runtime is sublinear as more kernels are added because the compiler can optimise each operation on a residue, rather than two for loops which would (probably?) lead to more cache misses.
  	It also looks very cool (in my opinion)!