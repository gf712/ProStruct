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

    int n_atoms = xyz.n_cols;

    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
#pragma omp simd reduction(+: sum_x, sum_y, sum_z)
    for (arma::uword i = 0; i < xyz.n_cols; ++i) {
        sum_x += xyz.at(0, i);
        sum_y += xyz.at(1, i);
        sum_z += xyz.at(2, i);
    }
    centroid[0] = sum_x / n_atoms;
    centroid[1] = sum_y / n_atoms;
    centroid[2] = sum_z / n_atoms;
}


void reposition_molecule(arma::mat &xyz, const arma::vec& centroid) {

    double x = centroid.at(0);
    double y = centroid.at(1);
    double z = centroid.at(2);

    for (arma::uword i = 0; i < xyz.n_cols; ++i) {
        xyz.at(0, i) -= x;
        xyz.at(1, i) -= y;
        xyz.at(2, i) -= z;
    }

}


void recentre_molecule(arma::mat &xyz) {

    arma::vec centroid(3);

    get_centroid(xyz, centroid);

    reposition_molecule(xyz, centroid);

}


void recentre_molecule(arma::mat &xyz, arma::mat &result) {

    arma::vec centroid(3);

    get_centroid(xyz, centroid);

    double x = centroid.at(0);
    double y = centroid.at(1);
    double z = centroid.at(2);

    for (arma::uword i = 0; i < xyz.n_cols; ++i) {
        result.at(0, i) -= x;
        result.at(1, i) -= y;
        result.at(2, i) -= z;
    }
}


void kabsch_rotation(arma::mat &xyz, arma::mat &other_xyz) {

    // source: https://en.wikipedia.org/wiki/Kabsch_algorithm

    arma::mat xyz_1(3, xyz.n_rows);
    arma::mat xyz_2(3, xyz.n_rows);

    recentre_molecule(xyz, xyz_1);
    recentre_molecule(other_xyz, xyz_2);

    arma::mat A = xyz_1.t() * xyz_2;

    arma::mat U;
    arma::vec s;
    arma::mat V;

    // First, calculate the SVD of the covariance matrix A.
    arma::svd(U, s, V, A);

    // Next, decide whether we need to correct our rotation matrix to ensure a right-handed coordinate system
    double d = arma::det(V * U.t()) > 0 ? 1 : -1;

    arma::mat I = arma::eye(3, 3);

    I.at(2, 2) = d;

    // Finally, calculate our optimal rotation matrix -> rotationMatrix
    arma::mat rotationMatrix = V * I * U.t();

    arma::mat P = xyz * U;

    arma::vec centroid_other(3);

    get_centroid(other_xyz, centroid_other);

    // centre xyz on other_xyz
    reposition_molecule(xyz, centroid_other);
}