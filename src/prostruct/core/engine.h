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
#include <prostruct/pdb/PDB.h>
#include <prostruct/utils/type_traits.h>

namespace prostruct::core {
	/**
	 * A near zero cost abstraction engine to execute multiple lambdas
	 * per residue in a loop.
	 *
	 */
	template <typename T, typename... Args>
	arma::Mat<T> residue_kernel_engine(
		const prostruct::residueVector<T>& residues, std::size_t start,
		Args... computations)
	{
		constexpr int n_computations = sizeof...(computations);
		auto result
			= arma::Mat<T>(n_computations, residues.size(), arma::fill::zeros);

		std::tuple<Args...> comp{ computations... };

		if constexpr (utils::lambdas_have_same_arity<Args...>()) {
			// the decay_t is necessary to remove information about where the
			// lambda was defined i.e. if defined in main, std::get<0>(comp)
			// would have type main()::lamda(T arg, ...) the decay_t returns
			// lambda(T arg, ...)
			constexpr size_t window_size = utils::lambda_properties<
				std::decay_t<decltype(std::get<0>(comp))>>::size;
#pragma omp parallel for
			for (std::size_t i = start; i < residues.size() - window_size + 1;
				 ++i) {
				execute_tuple(comp,
					vector_to_tuple_helper(residues,
						std::make_index_sequence<window_size>{}, i - start),
					result.col(i - start));
			}
		} else {
			static_assert(utils::lambdas_have_same_arity(computations...),
				"Not implemented yet!");
		}

		return result;
	}
}

#endif // PROSTRUCT_ENGINE_H
