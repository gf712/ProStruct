/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include <prostruct/core/kernels.h>
#include <prostruct/pdb/geometry.h>

namespace prostruct
{
	namespace geometry
	{

		template <typename T>
		T rmsd(const arma::Mat<T>& xyz, const arma::Mat<T>& xyz_other)
		{

			T sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
			for (arma::uword i = 0; i < xyz.n_cols; ++i)
			{
				sum += kernels::distance_lazy<T>(xyz.col(i), xyz_other.col(i));
			}

			return std::sqrt(sum / xyz.n_cols);
		}

		template <typename T>
		void get_centroid(const arma::Mat<T>& xyz, arma::Col<T>& centroid)
		{
			// the centroid is just the average point in cartesian point
			// does not take into account anything but coordinates

			size_t n_atoms = xyz.n_cols;

			T sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
#pragma omp simd reduction(+ : sum_x, sum_y, sum_z)
			for (arma::uword i = 0; i < xyz.n_cols; ++i)
			{
				sum_x += xyz.at(0, i);
				sum_y += xyz.at(1, i);
				sum_z += xyz.at(2, i);
			}
			centroid[0] = sum_x / n_atoms;
			centroid[1] = sum_y / n_atoms;
			centroid[2] = sum_z / n_atoms;
		}

		template <typename T>
		void reposition_molecule(arma::Mat<T>& xyz, const arma::Col<T>& centroid)
		{

			T x = centroid.at(0);
			T y = centroid.at(1);
			T z = centroid.at(2);

			for (arma::uword i = 0; i < xyz.n_cols; ++i)
			{
				xyz.at(0, i) -= x;
				xyz.at(1, i) -= y;
				xyz.at(2, i) -= z;
			}
		}

		template <typename T>
		void recentre_molecule(arma::Mat<T>& xyz)
		{

			arma::Col<T> centroid(3);

			get_centroid(xyz, centroid);

			reposition_molecule(xyz, centroid);
		}

		template <typename T>
		void recentre_molecule(arma::Mat<T>& xyz, arma::Mat<T>& result)
		{

			arma::Col<T> centroid(3);

			get_centroid(xyz, centroid);

			// copy array
			result = xyz;

			reposition_molecule(result, centroid);
		}

		template <typename T>
		void kabsch_rotation_(arma::Mat<T>& xyz, arma::Mat<T>& other_xyz)
		{

			// source: https://en.wikipedia.org/wiki/Kabsch_algorithm

			arma::Mat<T> xyz_1(3, xyz.n_cols);
			arma::Mat<T> xyz_2(3, xyz.n_cols);

			recentre_molecule(xyz, xyz_1);
			recentre_molecule(other_xyz, xyz_2);

			arma::Mat<T> A = xyz_1 * xyz_2.t();

			arma::Mat<T> U;
			arma::Col<T> s;
			arma::Mat<T> V;

			// First, calculate the SVD of the covariance matrix A.
			arma::svd(U, s, V, A);

			// Next, decide whether we need to correct our rotation matrix to
			// ensure a right-handed coordinate system
			arma::Mat<T> I = arma::eye<arma::Mat<T>>(3, 3);
			I.at(2, 2) = arma::det(V * U.t()) > 0 ? 1 : -1;

			// Finally, calculate our optimal rotation matrix -> rotationMatrix
			arma::Mat<T> rotationMatrix = V * I * U.t();

			// and apply rotation to the coordinate system
			xyz = rotationMatrix * xyz;
		}

		template <typename T>
		T kabsch_rmsd_(arma::Mat<T>& xyz, arma::Mat<T>& other_xyz)
		{

			kabsch_rotation_(xyz, other_xyz);

			return rmsd(xyz, other_xyz);
		}

		template float rmsd(const arma::Mat<float>&, const arma::Mat<float>&);

		template double rmsd(const arma::Mat<double>&, const arma::Mat<double>&);

		template void get_centroid(const arma::Mat<float>&, arma::Col<float>&);

		template void get_centroid(const arma::Mat<double>&, arma::Col<double>&);

		template float kabsch_rmsd_(arma::Mat<float>&, arma::Mat<float>&);
		template double kabsch_rmsd_(arma::Mat<double>&, arma::Mat<double>&);

		template void kabsch_rotation_(arma::Mat<float>&, arma::Mat<float>&);
		template void kabsch_rotation_(arma::Mat<double>&, arma::Mat<double>&);

		template void recentre_molecule(arma::Mat<float>&);

		template void recentre_molecule(arma::Mat<double>&);

	}
}