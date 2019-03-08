/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include <prostruct/struct/chain.h>

#include <iostream>

using namespace prostruct;

template <typename T>
Chain<T>::Chain(
	std::vector<std::shared_ptr<Residue<T>>> residues, std::string chainName)
{

	bool first = true;
	this->m_residues = residues;
	std::shared_ptr<Residue<T>> previousResidue;
	this->m_natoms = 0;

	// link residues
	for (auto const& residue : residues) {

		//        std::cout << residue->getResidueName() << std::endl;

		if (first) {
			previousResidue = residue;
			first = false;
		}

		else {
			residue->link(previousResidue);
			previousResidue = residue;
		}

		this->m_natoms += residue->n_atoms();
	}

	this->m_nresidues = static_cast<int>(residues.size());
}

template class prostruct::Chain<float>;
template class prostruct::Chain<double>;