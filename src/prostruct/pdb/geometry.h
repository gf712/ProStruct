/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include "prostruct/struct/residue.h"
#include "prostruct/utils/tuple_utils.h"
#include <armadillo>
#include <array>

#ifndef DNDEBUG
#define ARMA_NO_DEBUG
#endif

namespace prostruct
{
	namespace geometry
	{
		template <typename T>
		void dssp(
			const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&);

		template <typename T>
		void kabsch_sander(const arma::Mat<T>&, arma::Mat<T>&);

		template <typename T>
		void predict_H_coords(arma::Mat<T>& H_coords, const arma::Mat<T>& C_coords,
			const arma::Mat<T>& O_coords, const arma::Mat<T>& N_coords);

		template <typename T>
		void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii, arma::Col<T>& asa,
			arma::uword n_atoms, T probe, arma::uword n_sphere_points);

		template <typename T>
		void get_neighbours(const arma::Mat<T>&, arma::Mat<T>&, int, const arma::Col<T>&);

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
		void recentre_molecule(const arma::Mat<T>& xyz, arma::Mat<T>& result);

		template <typename T>
		T dihedrals(const arma::Col<T>& atom1, const arma::Col<T>& atom2, const arma::Col<T>& atom3,
			const arma::Col<T>& atom4, T coef);
	}
}
#endif // PROSTRUCT_GEOMETRY_H
