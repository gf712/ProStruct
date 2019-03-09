/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include <prostruct/pdb/struct_base.h>
#include <prostruct/struct/residue.h>
#include <prostruct/utils/io.h>

#include <fmt/format.h>

namespace prostruct
{
	template <typename T>
	class PDB;

	template <typename T>
	class StructBase;

	template <typename T>
	class Chain : public StructBase<T>
	{

	public:
		Chain(const std::vector<std::shared_ptr<Residue<T>>>&, const std::string&);
		Chain(const std::vector<std::shared_ptr<Residue<T>>>&, const std::string&,
			const arma::Mat<T>& xyz);

		virtual std::string to_string()
		{
			return format(fmt("<prostruct.Chain {} precision, with {} atoms, {} "
							  "residues at {}>"),
				demangled_type<T>(), this->m_natoms, this->m_nresidues, fmt::ptr(this));
		}

#ifndef SWIG
		template <typename... Args>
		arma::Col<arma::uword> get_atom_indices(Args... patterns)
		{
			arma::Col<arma::uword> max_result(this->m_nresidues * this->m_natoms);
			arma::uword i = 0;
			arma::uword position_in_chain = 0;

			for (const auto& residue : this->m_residues)
			{
				arma::Col<arma::uword> residue_atoms
					= residue->get_atom_indices(std::forward<Args>(patterns)...);
				max_result.subvec(arma::span(i, i + residue_atoms.n_rows - 1))
					= residue_atoms + position_in_chain;
				i += residue_atoms.n_rows;
				position_in_chain += residue->n_atoms();
			}

			arma::Col<arma::uword> result = max_result.head(i);

			return result;
		}
#endif
	private:
		std::string m_chain_name;
	};
}

#endif // PROSTRUCT_CHAIN_H
