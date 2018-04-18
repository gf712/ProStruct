//
// Created by gil on 10/04/18.
//

#include "../include/geometry.h"

double rmsd(const arma::mat& xyz, const arma::mat& xyz_other) {

    double sum = 0.0;

#pragma omp parallel for reduction(+: sum)
    for (arma::uword i = 0; i < xyz.n_cols; ++i) {
        sum += arma::norm(xyz.col(i) - xyz_other.col(i), 2);
    }

    return std::sqrt(sum / xyz.n_cols);
}
