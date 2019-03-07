//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include "prostruct/struct/residue.h"
#include "prostruct/utils/tuple_utils.h"
#include <armadillo>
#include <array>

#ifndef DNDEBUG
#define ARMA_NO_DEBUG
#endif

namespace prostruct {
	namespace geometry {
		template <typename T>
		void dssp(const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&,
			const arma::Mat<T>&);

		template <typename T>
		void kabsch_sander(const arma::Mat<T>&, arma::Mat<T>&);

		template <typename T>
		void predict_H_coords(arma::Mat<T>& H_coords,
			const arma::Mat<T>& C_coords, const arma::Mat<T>& O_coords,
			const arma::Mat<T>& N_coords);

		template <typename T>
		void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii,
			arma::Col<T>& asa, arma::uword n_atoms, T probe,
			arma::uword n_sphere_points);

		template <typename T>
		void get_neighbours(
			const arma::Mat<T>&, arma::Mat<T>&, int, const arma::Col<T>&);

		template <typename T>
		T rmsd(const arma::Mat<T>& xyz, const arma::Mat<T>& xyz_other);

		template <typename T>
		T kabsch_rmsd_(arma::Mat<T>& xyz, arma::Mat<T>& xyz_other);

		template <typename T>
		void kabsch_rotation_(arma::Mat<T>& xyz, arma::Mat<T>& xyz_other);
		template <typename T>
		void get_centroid(const arma::Mat<T>& xyz, arma::Col<T>& centroid);

		template <typename T> void recentre_molecule(arma::Mat<T>& xyz);

		template <typename T>
		void recentre_molecule(const arma::Mat<T>& xyz, arma::Mat<T>& result);

		template <typename T>
		T dihedrals(const arma::Col<T>& atom1, const arma::Col<T>& atom2,
			const arma::Col<T>& atom3, const arma::Col<T>& atom4, T coef);

		/**
		 *	A near zero cost abstraction engine to execute multiple lambdas
		 * per residue in a loop.
		 *
		 */
		template <size_t window_size, typename T, typename... Args>
		arma::Mat<T> atom_calculation_engine(
			const prostruct::residueVector<T>& residues, std::size_t start,
			Args... computations)
		{
			constexpr int n_computations = sizeof...(computations);
			auto result = arma::Mat<T>(
				n_computations, residues.size(), arma::fill::zeros);

			std::tuple<Args...> comp { computations... };
#pragma omp parallel for
			for (std::size_t i = start; i < residues.size() - window_size + 1;
				 ++i) {
				execute_tuple(comp,
					vector_to_tuple_helper(residues,
						std::make_index_sequence<window_size> {}, i - start),
					result.col(i - start));
			}

			return result;
		}

		/**
		 * An untyped version of dihedrals for armadillo optimisations
		 */
		template <typename T1, typename T2>
		inline T2 dihedrals_lazy(
			T1&& atom1, T1&& atom2, T1&& atom3, T1&& atom4, T2 coef)
		{
			arma::Col<T2> b1 = arma::normalise(atom1 - atom2);
			arma::Col<T2> b2 = arma::normalise(atom2 - atom3);
			arma::Col<T2> b3 = arma::normalise(atom3 - atom4);
			arma::Col<T2> n1 = arma::cross(b1, b2);
			arma::Col<T2> n2 = arma::cross(b2, b3);
			return std::atan2(
					   arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2))
				* coef;
		}
	}
}
#endif // PROSTRUCT_GEOMETRY_H
