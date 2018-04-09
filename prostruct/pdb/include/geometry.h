//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

void dssp(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&);
void kabsch_sander(const arma::mat&, const arma::mat&, const arma::mat&, const arma::mat&,
                   std::vector<bool>&, arma::mat&, const arma::uword);


#endif //PROSTRUCT_GEOMETRY_H
