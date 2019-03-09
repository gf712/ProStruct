/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "prostruct/pdb/geometry.h"
#include "prostruct/struct/atom.h"

namespace prostruct
{
	namespace geometry
	{

		constexpr double golden_angle = 2.399963229728653;

		template <typename T>
		void generate_sphere(arma::Mat<T>& result)
		{

			T offset = 2.0 / result.n_cols;

			for (arma::uword i = 0; i < result.n_cols; ++i)
			{

				T y = i * offset - 1.0 + (offset / 2.0);
				T r = std::sqrt(1.0 - y * y);

				T t = i * golden_angle;

				result.at(0, i) = r * std::cos(t); // x
				result.at(1, i) = y; // y
				result.at(2, i) = r * std::sin(t); // z
			}
		}

		template <typename T>
		void calculate_atom_SASA(const arma::Mat<T>& xyz, const arma::Col<T>& radius,
			const arma::subview_col<T>& neighbours, const arma::uword current_atom_index,
			const T probe, const arma::Mat<T>& sphere_points, const T adjustment, arma::Col<T>& asa)
		{

			arma::Col<T> atom_XYZ = xyz.col(current_atom_index);
			arma::uvec neighbourIndices = arma::find(neighbours, 0);

			T atomRadius = probe + radius.at(current_atom_index);
			arma::uword nNeighbours = neighbourIndices.size();
			arma::uword accessiblePoints = 0;

			arma::uword k = 0;

			for (arma::uword i = 0; i < sphere_points.n_cols; ++i)
			{
				for (arma::uword j = k; j < nNeighbours + k; ++j)
				{

					arma::uword index = neighbourIndices.at(j % nNeighbours);
					T r_2 = std::pow(radius.at(index) + probe, 2);
					T dist = arma::sum(
						arma::square((sphere_points(arma::span::all, i) * atomRadius + atom_XYZ)
							- xyz.col(index)));

					if (dist < r_2)
					{
						k = j;
						goto BURRIED;
					}
				}
				accessiblePoints++;
			BURRIED:
				continue;
			}

			asa.at(current_atom_index) = adjustment * accessiblePoints * std::pow(atomRadius, 2);
		}

		template <typename T>
		void get_neighbours(const arma::Mat<T>& xyz, arma::Mat<T>& neighbours, int n_atoms,
			const arma::Col<T>& radii)
		{
#pragma omp parallel for
			for (arma::uword i = 0; i < n_atoms; ++i)
			{
				for (arma::uword j = 1; j < n_atoms; ++j)
				{
					T cutoff = std::pow(radii.at(i) + radii.at(j), 2);
					if (arma::sum(arma::square(xyz(arma::span::all, i) - xyz(arma::span::all, j)))
						< cutoff)
					{
						neighbours.at(i, j) = 1;
						neighbours.at(j, i) = 1;
					}
				}
			}
		}

		template <typename T>
		void shrake_rupley(const arma::Mat<T>& xyz, const arma::Col<T>& radii, arma::Col<T>& asa,
			arma::uword n_atoms, T probe, arma::uword n_sphere_points)
		{

			arma::Mat<T> sphere_points(3, n_sphere_points);
			arma::Mat<T> neighbours(n_atoms, n_atoms);

			generate_sphere(sphere_points);

			get_neighbours(xyz, neighbours, n_atoms, radii);

			T adjustment = 4.0 * M_PI / n_sphere_points;

#pragma omp parallel for
			for (arma::uword i = 0; i < n_atoms; ++i)
			{
				calculate_atom_SASA(
					xyz, radii, neighbours.col(i), i, probe, sphere_points, adjustment, asa);
			}
		}

		template void shrake_rupley(const arma::Mat<float>&, const arma::Col<float>&,
			arma::Col<float>&, arma::uword n_atoms, float probe, arma::uword n_sphere_points);

		template void shrake_rupley(const arma::Mat<double>&, const arma::Col<double>&,
			arma::Col<double>&, arma::uword n_atoms, double probe, arma::uword n_sphere_points);

	}
}