//
// Created by gil on 09/04/18.
//

#include "../include/geometry.h"

void dihedrals(const arma::cube &atoms, arma::vec &angles) {

    // atom coordinates are stored in a cube (3D tensor)
    // atoms has shape (3, 4, n_angles)
    // where n_angles is the number of angles to calculate

    for (arma::uword i = 0; i < atoms.n_slices; i++) {
        // access cube using slice -> returns matrix
        const arma::mat &atoms_i = atoms.slice(i);

        arma::vec vec1 = atoms_i.col(0) - atoms_i.col(1);
        arma::vec vec2 = atoms_i.col(2) - atoms_i.col(3);

        double num = arma::dot(vec1, vec2);
        double denom = arma::norm(vec1, 2) * arma::norm(vec2, 2);
        angles.at(i) = std::acos(num /denom);
        i++;
    }

}