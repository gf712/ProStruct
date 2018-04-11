//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include <armadillo>

void dssp(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&);
void kabsch_sander(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&,
                   std::vector<bool>&, arma::mat&, const arma::uword);
void shrake_rupley(const arma::mat &xyz, const arma::vec& radii, arma::vec &asa, int n_atoms, double probe);
double rmsd(const arma::mat& xyz, const arma::mat& xyz_other);


#endif //PROSTRUCT_GEOMETRY_H
