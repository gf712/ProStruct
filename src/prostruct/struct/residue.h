/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_RESIDUE_H
#define PROSTRUCT_RESIDUE_H

#include <prostruct/struct/atom.h>
#include <prostruct/struct/utils.h>
#include <prostruct/utils/type_traits.h>

namespace prostruct
{

	using aminoAcidAtomMap = std::vector<std::map<std::string, std::string>>;
	using stringIndexMap = std::map<std::string, int>;

	template <typename T>
	class Residue
	{
		template <typename T1>
		friend class Chain;

	public:
		// takes an arbitrary number of atoms and tries to form a residue
		//    Residue(std::unique_ptr<Atom> atoms...);
		Residue(atomVector<T>, const std::string&, const std::string&, bool = false, bool = false);

		inline static const std::vector<std::string> backbone_atom_names = { "C", "CA", "N", "O" };

		atomVector<T> getBackbone() const noexcept
		{
			std::vector<std::shared_ptr<Atom<T>>> backboneAtoms;
			for (auto const& i : backbone)
			{
				backboneAtoms.emplace_back(atoms[i]);
			}
			return backboneAtoms;
		}

		AminoAcid get_amino_acid_type() const noexcept { return m_amino_acid; }

		arma::Mat<T> get_backbone_atoms() const noexcept
		{
			return xyz(arma::span::all, arma::span(0, 3));
		}

		atomVector<T> get_sidechain() const noexcept
		{
			std::vector<std::shared_ptr<Atom<T>>> sidechainAtoms;
			for (auto const& i : sidechain)
			{
				sidechainAtoms.emplace_back(atoms[i]);
			}
			return sidechainAtoms;
		}

		arma::Mat<T> get_sidechain_atoms() const noexcept
		{
			if (get_amino_acid_type() == AminoAcid::GLY)
				return arma::Mat<T>(3, 0);
			else
				return xyz(arma::span::all, arma::span(4, atoms.size() - 1));
		}

		std::shared_ptr<Atom<T>> operator[](const int index) const { return atoms[index]; }

		arma::Mat<T> get_xyz() const noexcept { return xyz; }

		std::string get_name() const noexcept { return m_residue_name; }

		std::shared_ptr<Atom<T>> get_atom(int index) const noexcept { return atoms[index]; }

		void link(std::shared_ptr<Residue<T>>);

		void createBonds();

		int n_atoms() const noexcept { return atoms.size(); };

		atomVector<T> getAtoms() const noexcept { return atoms; }

		arma::Col<T> getRadii() const noexcept { return radii; }

		bool is_n_terminus() const noexcept { return m_n_terminus; }

		bool is_c_terminus() const noexcept { return m_c_terminus; }

#ifndef SWIG
		template <typename expect_one = std::true_type, typename... Args>
		arma::Col<arma::uword> get_atom_indices(const Args&... patterns) const noexcept
		{
			arma::Col<arma::uword> max_result;
			if constexpr (expect_one::value)
				max_result.set_size(sizeof...(Args));
			else
				max_result.set_size(atoms.size());

			std::vector<std::string> vec = { patterns... };
			arma::uword result_idx = 0;

			for (auto atom = atoms.cbegin(); atom != atoms.cend(); ++atom)
			{
				const auto pattern_result = std::find_if(vec.cbegin(), vec.cend(),
					[&atom](const auto& pattern) { return pattern == atom->get()->get_name(); });

				if (pattern_result == vec.cend())
					continue;
				else
				{
					max_result(result_idx) = std::distance(atoms.cbegin(), atom);
					++result_idx;
					if constexpr (expect_one::value)
					{
						if (result_idx == sizeof...(Args))
							break;
					}
				}
			}
			if constexpr (expect_one::value)
				return max_result;
			else
				return max_result.head(result_idx);
		}

		template <typename expect_one = std::true_type, typename... Args>
		arma::Mat<T> get_atom_coords(Args... idx) const noexcept
		{
			if constexpr (utils::all_integral_v<Args...>)
			{
				std::vector<arma::uword> idx_vector = { idx... };
				return xyz(arma::span::all, idx_vector);
			}
			if constexpr (utils::all_same_v<Args...>)
			{
				arma::Mat<T> max_result;
				if constexpr (expect_one::value)
					max_result.set_size(3, sizeof...(Args));
				else
					max_result.set_size(3, atoms.size());

				std::vector<std::string> vec = { idx... };
				arma::uword result_idx = 0;

				for (auto atom = atoms.cbegin(); atom != atoms.cend(); ++atom)
				{
					const auto pattern_result
						= std::find_if(vec.cbegin(), vec.cend(), [&atom](const auto& pattern) {
							  return pattern == atom->get()->get_name();
						  });

					if (pattern_result == vec.cend())
						continue;
					else
					{
						max_result.col(result_idx) = xyz.col(std::distance(atoms.cbegin(), atom));
						++result_idx;
						if constexpr (expect_one::value)
						{
							if (result_idx == sizeof...(Args))
								break;
						}
					}
				}
				if constexpr (expect_one::value)
					return max_result;
				else
					return max_result.head(result_idx);
			}
		}

		arma::Mat<T> get_atom_coords(const arma::Col<arma::uword>& idx) const noexcept
		{
			return xyz(arma::Col<arma::uword>({ 0, 1, 2 }), idx);
		}
#endif

	private:
		arma::Mat<T> xyz;
		arma::Col<T> radii;
		bool m_n_terminus;
		bool m_c_terminus;
		std::vector<int> backbone; /**< A vector with the index number of the
									  backbone atoms */
		std::vector<int> sidechain; /**< A vector with the index number of the
									   sidechain atoms */
		std::string aminoAcidName; /**< Name of the amino acid, e.g. ALA */
		enum AminoAcid m_amino_acid; /**< Amino acid enum, e.g. ALA */
		std::string m_residue_name; /**< Name of the residue, e.g. ALA1 */
		atomVector<T> atoms; /**< A vector with the pointers to the Atom objects */
		std::map<std::string, int> atomMap; /**< Map atom name to internal index */
	};

	template <typename T>
	using residueVector = std::vector<std::shared_ptr<Residue<T>>>;
}

#endif // PROSTRUCT_RESIDUE_H
