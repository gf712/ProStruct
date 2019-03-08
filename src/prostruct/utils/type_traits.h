/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_TYPE_TRAITS_H
#define PROSTRUCT_TYPE_TRAITS_H

#include <tuple>
#include <type_traits>

namespace prostruct::utils {
	// adapted from https://qiita.com/angeart/items/94734d68999eca575881
	namespace detail {
		template <typename ReturnType, typename Cls, typename... Args>
		struct lambda_properties_helper {
			using return_type = ReturnType;
			static constexpr size_t size = sizeof...(Args) - 1;
			template <size_t i> struct arg {
				typedef
					typename std::tuple_element<i, std::tuple<Args...>>::type
						type;
			};
		};

	}

	template <typename LambdaType>
	struct lambda_properties : lambda_properties<decltype(&LambdaType::operator())>
	{
	};

	template <typename ReturnType, typename Cls, typename... Args>
	struct lambda_properties<ReturnType (Cls::*)(Args...)>
		: detail::lambda_properties_helper<ReturnType, Cls, std::true_type, Args...>
	{
	};

	template <typename ReturnType, typename Cls, typename... Args>
	struct lambda_properties<ReturnType (Cls::*)(Args...) const>
		: detail::lambda_properties_helper<ReturnType, Cls, std::false_type, Args...>
	{
	};

	template <typename Head, typename ...Tail>
	constexpr bool lambda_compare_arity()
	{
		return ((lambda_properties<Head>::size == lambda_properties<Tail>::size) && ...);
	}

	template <typename ...LambdaTypes>
	constexpr bool lambdas_have_same_arity()
	{
		if constexpr (sizeof...(LambdaTypes) == 0)
			return false;
		else if constexpr (sizeof...(LambdaTypes) == 1)
			return true;
		else
			return lambda_compare_arity<LambdaTypes...>();
	}

	template <typename ...Args>
	inline constexpr bool all_integral_v = (std::is_integral<Args>::value && ...);

	template <typename...>
	inline constexpr auto all_same_v = std::true_type{};

	template <typename T, typename... Rest>
	inline constexpr auto all_same_v<T, Rest...> = (std::is_same_v<T, Rest> && ...);
}

#endif // PROSTRUCT_TYPE_TRAITS_H
