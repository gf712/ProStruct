//
// Created by gil on 06/04/18.
//

#ifndef PROSTRUCT_GEOMETRY_H
#define PROSTRUCT_GEOMETRY_H

#include <armadillo>
#include <array>
#include <cmath>
#include "prostruct/struct/residue.h"
#include "prostruct/utils/tuple_utils.h"


#ifndef DNDEBUG
#define ARMA_NO_DEBUG
#endif

namespace prostruct {
	namespace geometry {
		template <typename T>
		void dssp(const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&, const arma::Mat<T>&);

		template <typename T>
		void kabsch_sander(const arma::Mat<T>&, arma::Mat<T>&);

		template <typename T>
		void predict_H_coords(arma::Mat<T>& H_coords, const arma::Mat<T>& C_coords, const arma::Mat<T>& O_coords,
			const arma::Mat<T>& N_coords);

		template <typename T>
		void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii, arma::Col<T>& asa, arma::uword n_atoms,
		        T probe, arma::uword n_sphere_points);

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
		void dihedrals(const arma::Cube<T>& atoms, arma::Col<T>& angles);

		/** 
		  *	A near zero cost abstraction engine to execute multiple lambdas
		  * per residue in a loop.
		  *
		  */
		template<size_t window_size, typename T, typename ...Args>
		arma::Mat<T> atom_calculation_engine(const prostruct::residueVector<T>& residues, Args... computations)
		{
			constexpr int n_computations = sizeof...(computations);
			auto result = arma::Mat<T>(n_computations, residues.size(), arma::fill::zeros);

			std::tuple<Args...> comp {computations...};

			int offset = std::ceil((window_size-1.0) / 2.0);
			int remainder = window_size - offset;

			for (int i = offset; i < residues.size() - remainder; ++i)
			{
				execute_tuple(comp, vector_to_tuple_helper(residues, std::make_index_sequence<window_size>{}, i - offset), result.col(i - offset));
			}

			return result;
		}
	}
}
#endif //PROSTRUCT_GEOMETRY_H
