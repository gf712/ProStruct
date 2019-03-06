//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include "prostruct/struct/residue.h"

template <typename T>
using chainAtomVector = std::vector<atomVector<T>>;

namespace prostruct {
	template <typename T>
	class Chain {
		template <typename T1>
	    friend class PDB;

	public:
		Chain(std::vector<std::shared_ptr<Residue<T>>>, std::string);
		int n_residues() { return m_nresidues; };
		int n_atoms() { return m_natoms; };
		residueVector<T> getResidues() { return residues; }
		chainAtomVector<T> getBackboneAtoms();

	private:
        std::string chainName;
		std::vector<std::shared_ptr<Residue<T>>> residues;
		int m_nresidues;
		int m_natoms;
	protected:
		template <typename... Args>
		arma::Col<arma::uword> get_atom_indices(Args... patterns)
		{
			arma::Col<arma::uword> max_result(m_nresidues*m_natoms);
			arma::uword i = 0;
			arma::uword position_in_chain = 0;

		    for(const auto& residue: residues) {
				arma::Col<arma::uword> residue_atoms = residue->get_atom_indices(std::forward<Args>(patterns)...);
		        max_result.subvec(arma::span(i, i+residue_atoms.n_rows-1)) = residue_atoms + position_in_chain;
		    	i+=residue_atoms.n_rows;
		        position_in_chain += residue->n_atoms();
		    }

		    arma::Col<arma::uword> result = max_result.head(i);

			return result;
		}
	};
}

#endif //PROSTRUCT_CHAIN_H
