//
// Created by gil on 27/03/18.
//

#include <iostream>
#include "prostruct/struct/chain.h"

template <typename T>
Chain<T>::Chain(std::vector<std::shared_ptr<Residue<T>>> residues_, std::string chainName_) {

    bool first = true;
    residues = residues_;
    std::shared_ptr<Residue<T>> previousResidue;
    nAtoms = 0;

    // link residues
    for (auto const &residue: residues_) {

//        std::cout << residue->getResidueName() << std::endl;

        if (first) {
            previousResidue = residue;
            first = false;
        }

        else {
            residue->link(previousResidue);
            previousResidue = residue;
        }

        nAtoms += residue->n_atoms();
    }

    nResidues = static_cast<int>(residues.size());
}

template <typename T>
chainAtomVector<T> Chain<T>::getBackboneAtoms() {
    std::vector<std::vector<std::shared_ptr<Atom<T>>>> backboneAtoms;
    backboneAtoms.reserve((nResidues));
    for(auto const& residue: residues) {
        backboneAtoms.push_back(residue->getBackbone());
    }
    return backboneAtoms;
}

template class Chain<float>;
template class Chain<double>;