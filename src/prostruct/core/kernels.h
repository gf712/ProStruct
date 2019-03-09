/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_KERNELS_H
#define PROSTRUCT_KERNELS_H

#include <armadillo>

namespace prostruct::kernels
{
	/**
	 * An untyped version of dihedrals for armadillo optimisations
	 */
	template <typename T1, typename T2>
	inline T2 dihedrals_lazy(T1&& atom1, T1&& atom2, T1&& atom3, T1&& atom4, T2 coef)
	{
		arma::Col<T2> b1 = arma::normalise(atom1 - atom2);
		arma::Col<T2> b2 = arma::normalise(atom2 - atom3);
		arma::Col<T2> b3 = arma::normalise(atom3 - atom4);
		arma::Col<T2> n1 = arma::cross(b1, b2);
		arma::Col<T2> n2 = arma::cross(b2, b3);
		return std::atan2(arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2)) * coef;
	}
	/**
	 * An untyped version of the norm for armadillo
	 * optimisations.
	 */
	template <typename T1, typename T2>
	inline T1 norm_lazy(T2&& atom1, T2&& atom2)
	{
		return arma::dot(atom1 - atom2, atom1 - atom2);
	}

	/**
	 * An untyped version of euclidean distances for armadillo
	 * optimisations.
	 * Note that the dot product tends to be quicker
	 * than arma::norm when using the MKL backend.
	 */
	template <typename T1, typename T2>
	inline T1 distance_lazy(T2&& atom1, T2&& atom2)
	{
		return std::sqrt(norm_lazy<T1>(atom1, atom2));
	}
}

#endif // PROSTRUCT_KERNELS_H
