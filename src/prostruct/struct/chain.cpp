//
// Created by gil on 27/03/18.
//

#include "prostruct/struct/chain.h"
#include "chain.h"

#include <iostream>

using namespace prostruct;

template <typename T>
Chain<T>::Chain(
	std::vector<std::shared_ptr<Residue<T>>> residues_, std::string chainName_)
{

	bool first = true;
	residues = residues_;
	std::shared_ptr<Residue<T>> previousResidue;
	m_natoms = 0;

	// link residues
	for (auto const& residue : residues_) {

		//        std::cout << residue->getResidueName() << std::endl;

		if (first) {
			previousResidue = residue;
			first = false;
		}

		else {
			residue->link(previousResidue);
			previousResidue = residue;
		}

		m_natoms += residue->n_atoms();
	}

	m_nresidues = static_cast<int>(residues.size());
}

template <typename T> chainAtomVector<T> Chain<T>::getBackboneAtoms()
{
	std::vector<std::vector<std::shared_ptr<Atom<T>>>> backboneAtoms;
	backboneAtoms.reserve(m_nresidues);
	for (auto const& residue : residues) {
		backboneAtoms.push_back(residue->getBackbone());
	}
	return backboneAtoms;
}

template class prostruct::Chain<float>;
template class prostruct::Chain<double>;