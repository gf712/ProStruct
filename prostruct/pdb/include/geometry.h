//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include <armadillo>

void dssp(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&);
void kabsch_sander(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&,
                   std::vector<bool>&, arma::mat&, const arma::uword);
static void predict_H_coords(arma::mat &H_coords, const arma::mat &C_coords, const arma::mat &O_coords, const arma::mat &N_coords);
void shrake_rupley(const arma::mat &xyz, const arma::vec& radii, arma::vec &asa, int n_atoms, double probe);
double rmsd(const arma::mat& xyz, const arma::mat& xyz_other);
void get_centroid(const arma::mat &xyz, arma::vec &centroid);
void recentre_molecule(arma::mat &xyz);
void recentre_molecule(arma::mat &xyz, arma::mat &result);
void dihedrals(const arma::cube &atoms, arma::vec &angles);

#endif //PROSTRUCT_GEOMETRY_H
