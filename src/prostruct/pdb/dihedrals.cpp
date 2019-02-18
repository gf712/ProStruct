//
// Created by gil on 09/04/18.
//

#include "prostruct/pdb/geometry.h"

template <typename T>
void dihedrals(const arma::Cube<T>& atoms, arma::Col<T>& angles)
{

	// atom coordinates are stored in a cube (3D tensor)
	// atoms has shape (3, 4, n_angles)
	// where n_angles is the number of angles to calculate

	for (arma::uword i = 0; i < atoms.n_slices; i++) {
		// access cube using slice -> returns matrix
		const arma::Mat<T>& atoms_i = atoms.slice(i);

		arma::Col<T> b1 = arma::normalise(atoms_i.col(0) - atoms_i.col(1));
		arma::Col<T> b2 = arma::normalise(atoms_i.col(1) - atoms_i.col(2));
		arma::Col<T> b3 = arma::normalise(atoms_i.col(2) - atoms_i.col(3));

		//        std::cout << "Atom: " << i << std::endl;
		//        std::cout << "b1" << std::endl;
		//        b1.print();
		//        std::cout << "b2" << std::endl;
		//        b2.print();
		//        std::cout << "b3" << std::endl;
		//        b3.print();

		arma::Col<T> n1 = arma::cross(b1, b2);
		arma::Col<T> n2 = arma::cross(b2, b3);

		//        std::cout << "n1" << std::endl;
		//        n1.print();
		//
		//        std::cout << "n2" << std::endl;
		//        n2.print();

		//        std::cout << "x" << std::endl;
		//        n1.print();
		//
		//        std::cout << "y" << std::endl;
		//        n2.print();

		angles.at(i) = std::atan2(arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2));

		//        std::cout << angles.at(i) * (180.0 / M_PI) << std::endl;
	}
}

template void dihedrals(const arma::Cube<float>&, arma::Col<float>&);
template void dihedrals(const arma::Cube<double>&, arma::Col<double>&);