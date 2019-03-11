/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_UTILS_H
#define PROSTRUCT_UTILS_H

#include <memory>
#include <vector>

namespace prostruct
{

	template <typename T>
	class Atom;

	template <typename T>
	class Residue;

	enum class aaLocation
	{
		Backbone,
		Sidechain
	};

	enum class AminoAcid
	{
		ARG,
		ALA,
		ASN,
		ASP,
		CYS,
		GLN,
		GLU,
		GLY,
		HIS,
		ILE,
		LEU,
		LYS,
		MET,
		PHE,
		PRO,
		SER,
		THR,
		TRP,
		TYR,
		VAL
	};

	template <typename T>
	using atomVector = std::vector<std::shared_ptr<Atom<T>>>;

	template <typename T>
	using residueVector = std::vector<std::shared_ptr<Residue<T>>>;
};

#endif // PROSTRUCT_UTILS_H
