/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "PDB.h"
#include <prostruct/pdb/PDB.h>

using namespace prostruct;

template <typename T>
PDB<T>::PDB(const std::string& filename)
	: StructBase<T>()
	, m_filename(filename)
{
	std::map<std::string, std::map<std::string, atomVector<T>, AASequenceOrder>> chainAtomMap;

	createMap(m_filename, chainAtomMap, m_chain_order);

	this->m_natoms = 0;
	this->m_nresidues = 0;
	this->m_xyz.set_size(3, 0);
	arma::uword start_current_atom = 0;
	arma::uword end_current_atom = 0;
	for (const auto& chain : m_chain_order)
	{
		residueVector<T> residues;
		const std::map<std::string, atomVector<T>, AASequenceOrder> chain_i
			= chainAtomMap.at(chain);
		this->append_new_residue(*chain_i.cbegin(), residues, true, false);
		for (auto atomPair = std::next(chain_i.cbegin()); atomPair != std::prev(chain_i.cend());
			 ++atomPair)
		{
			this->append_new_residue(*atomPair, residues, false, false);
		}
		this->append_new_residue(*std::prev(chain_i.cend()), residues, false, true);
		end_current_atom += this->m_natoms - start_current_atom;
		m_chain_map[chain] = std::make_shared<Chain<T>>(residues, chain,
			this->m_xyz(arma::span::all, arma::span(start_current_atom, end_current_atom - 1)));
		start_current_atom = end_current_atom;
		this->m_nresidues += static_cast<arma::uword>(residues.size());
		this->m_residues.insert(this->m_residues.end(), residues.begin(), residues.end());
	}
	this->m_number_of_chains = static_cast<int>(m_chain_map.size());
}

template <typename T>
PDB<T> PDB<T>::fetch(std::string PDB_id)
{
	throw "Not implemented";
}

// template <typename T>
// arma::Mat<T> PDB<T>::get_backbone_atoms() const noexcept
// {
// 	arma::Mat<T> result(3, this->m_nresidues * 4);
// 	arma::uword pos = 0;
// 	arma::uword i = 0;

// 	for (auto const& chain : m_chain_order)
// 	{
// 		for (auto const& residue : m_chain_map.at(chain)->get_residues())
// 		{
// 			result(arma::span::all, arma::span(i, i + 3))
// 				= this->m_xyz(arma::span::all, arma::span(pos, pos + 3));
// 			pos += residue->n_atoms();
// 			i += 4;
// 		}
// 	}
// 	return result;
// }

template class prostruct::PDB<float>;
template class prostruct::PDB<double>;