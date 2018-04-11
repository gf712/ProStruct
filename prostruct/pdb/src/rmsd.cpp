//
// Created by gil on 10/04/18.
//

#include "../include/geometry.h"

double rmsd(const arma::mat& xyz, const arma::mat& xyz_other) {

    return arma::sum(arma::norm(xyz - xyz_other, 2)) / xyz.n_cols;

}
