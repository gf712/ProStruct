/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_ENGINE_H
#define PROSTRUCT_ENGINE_H

#include <armadillo>

namespace prostruct::core {
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
		auto result
			= arma::Mat<T>(n_computations, residues.size(), arma::fill::zeros);

		std::tuple<Args...> comp { computations... };

		for (std::size_t i = start; i < residues.size() - window_size + 1;
			 ++i) {
			execute_tuple(comp,
				vector_to_tuple_helper(residues,
					std::make_index_sequence<window_size> {}, i - start),
				result.col(i - start));
		}

		return result;
	}
}

#endif // PROSTRUCT_ENGINE_H
