//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include <armadillo>

#ifndef DNDEBUG
#define ARMA_NO_DEBUG
#endif

template <typename T>
void dssp(const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&);

template <typename T>
void kabsch_sander(const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&,
	arma::Mat<T>&, const arma::uword);

template <typename T>
void predict_H_coords(arma::Mat<T>& H_coords, const arma::Mat<T>& C_coords, const arma::Mat<T>& O_coords, const arma::Mat<T>& N_coords);

template <typename T>
void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii, arma::Col<T>& asa, int n_atoms, T probe);

template <typename T>
T rmsd(const arma::Mat<T>& xyz, const arma::Mat<T>& xyz_other);

template <typename T>
T kabsch_rmsd_(arma::Mat<T>& xyz, arma::Mat<T>& xyz_other);

template <typename T>
void kabsch_rotation_(arma::Mat<T>& xyz, arma::Mat<T>& xyz_other);

template <typename T>
void get_centroid(const arma::Mat<T>& xyz, arma::Col<T>& centroid);

template <typename T>
void recentre_molecule(arma::Mat<T>& xyz);

template <typename T>
void recentre_molecule(arma::Mat<T>& xyz, arma::Mat<T>& result);

template <typename T>
void dihedrals(const arma::Cube<T>& atoms, arma::Col<T>& angles);

#endif //PROSTRUCT_GEOMETRY_H
