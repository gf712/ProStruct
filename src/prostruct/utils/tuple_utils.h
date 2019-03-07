/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include <tuple>
#include <type_traits>
#include <utility>
#include <armadillo>

namespace prostruct {
	// taken from https://stackoverflow.com/a/12650100
	template<size_t N>
	struct Apply {
	    template<typename F, typename T, typename... A>
	    static inline auto apply(F && f, T && t, A &&... a) {
	        return Apply<N-1>::apply(::std::forward<F>(f), ::std::forward<T>(t),
	            ::std::get<N-1>(::std::forward<T>(t)), ::std::forward<A>(a)...
	        );
	    }
	};

	template<>
	struct Apply<0> {
	    template<typename F, typename T, typename... A>
	    static inline auto apply(F && f, T &&, A &&... a) {
	        return ::std::forward<F>(f)(::std::forward<A>(a)...);
	    }
	};

	template<typename F, typename T, typename ResultType>
	void apply(F && f, T && t, ResultType& result) {
	    result = Apply< ::std::tuple_size< ::std::decay_t<T>
	      >::value>::apply(::std::forward<F>(f), ::std::forward<T>(t));
	}

	// adapted from https://en.cppreference.com/w/cpp/utility/integer_sequence
	template <typename LambdaTuple, typename ...LambdaArgs, std::size_t... Idx, typename ResultType>
	void execute_tuple_helper(const LambdaTuple& lambda_tuple, 
		const std::tuple<LambdaArgs...>& lambda_args, 
		std::index_sequence<Idx...>,
		ResultType&& result)
	{
		return (apply(std::get<Idx>(lambda_tuple), lambda_args, result[Idx]), ...);
	}

	template <typename ...Args, typename ...LambdaArgs, typename ResultType>
	void execute_tuple(const std::tuple<Args...>& lambda_tuple, 
		const std::tuple<LambdaArgs...>& lambda_args,
		ResultType&& result)
	{
		return execute_tuple_helper(lambda_tuple, lambda_args, std::index_sequence_for<Args...>{}, std::forward<ResultType>(result));
	}

	template <typename T, std::size_t... Idx>
	auto vector_to_tuple_helper(const std::vector<T>& vec, std::index_sequence<Idx...>, size_t offset) {
 		return std::make_tuple(vec[Idx + offset]...);
	}
}