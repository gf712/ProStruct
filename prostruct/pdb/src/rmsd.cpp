//
// Created by gil on 10/04/18.
//

#include "../include/geometry.h"

double rmsd(const arma::mat& xyz, const arma::mat& xyz_other) {

    double sum = 0.0;

    for (int i = 0; i < xyz.n_cols; ++i) {
        sum += arma::norm(xyz - xyz_other, 2);
    }

    return sum / xyz.n_cols;

}
