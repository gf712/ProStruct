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


void get_centroid(const arma::mat &xyz, arma::vec &centroid) {
    // the centroid is just the average point in cartesian point
    // does not take into account anything but coordinates

    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
#pragma omp simd reduction(+: sum_x, sum_y, sum_z)
    for (arma::uword i = 0; i < xyz.n_cols; ++i) {
        sum_x += xyz.at(0, i);
        sum_y += xyz.at(1, i);
        sum_z += xyz.at(2, i);
    }
    centroid[0] = sum_x;
    centroid[1] = sum_y;
    centroid[2] = sum_z;

    centroid /= xyz.n_cols;
}


void recentre_molecules(const arma::mat &xyz, const arma::mat &xyz_other) {



}
