//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include "residue.h"

template <typename T>
using chainAtomVector = std::vector<atomVector<T>>;

template <typename T>
class Chain {
public:
    Chain(std::vector<std::shared_ptr<Residue<T>>>, std::string);
    int n_residues() { return nResidues; };
    int n_atoms() { return nAtoms; };
    std::vector<std::shared_ptr<Residue<T>>> getResidues() { return residues; }
    chainAtomVector<T> getBackboneAtoms();

private:
    std::string chainName;
    std::vector<std::shared_ptr<Residue<T>>> residues;
    int nResidues;
    int nAtoms;
};


#endif //PROSTRUCT_CHAIN_H
