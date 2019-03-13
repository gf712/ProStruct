/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#include <prostruct/parsers/PDBparser.h>
#include <prostruct/pdb/struct_base.h>
#include <prostruct/struct/chain.h>
#include <prostruct/utils/io.h>
#include <prostruct/utils/ast.h>

#include <set>

#include <fmt/format.h>

namespace prostruct
{

	template <class T>
	class StructBase;
	template <typename T>
	class Chain;

	template <typename T>
	class PDB : public StructBase<T>
	{
	public:
		PDB(const std::string& filename);

		static PDB fetch(std::string);

		virtual std::string to_string() const
		{
			return format(fmt("<prostruct.PDB {} precision, with {} atoms, {} "
							  "residues at {}>"),
				demangled_type<T>(), this->m_natoms, this->m_nresidues, fmt::ptr(this));
		}

		std::string get_filename() { return m_filename; }

		std::shared_ptr<Chain<T>> get_chain(const std::string& name) const
		{
			return m_chain_map.at(name);
		}

		std::vector<std::string> get_chain_names() const noexcept { return m_chain_order; }

#ifndef SWIG
		template <typename... Args>
		arma::Col<arma::uword> get_atom_indices(Args... patterns) const noexcept
		{
			arma::Col<arma::uword> max_result(this->m_nresidues * this->m_natoms);
			arma::uword i = 0;
			arma::uword position_in_pdb = 0;

			for (const auto& chain_name : m_chain_order)
			{
				arma::Col<arma::uword> residue_atoms
					= m_chain_map[chain_name]->get_atom_indices(std::forward<Args>(patterns)...);
				max_result(arma::span(i, i + residue_atoms.n_rows - 1))
					= residue_atoms + position_in_pdb;
				i += residue_atoms.n_rows;
				position_in_pdb += m_chain_map[chain_name]->m_natoms;
			}

			arma::Col<arma::uword> result = max_result.head(i);

			return result;
		}
#endif
		// virtual arma::Mat<T> get_backbone_atoms() const noexcept override;

		int n_chains() { return m_number_of_chains; }

	private:
		std::string m_filename;
		std::vector<std::string> m_chain_order;
		std::map<std::string, std::shared_ptr<Chain<T>>> m_chain_map;
		int m_number_of_chains;
	};

} // namespace prostruct

#endif // PROSTRUCT_PDB_H
