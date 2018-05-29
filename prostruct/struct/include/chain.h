//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include "residue.h"

typedef std::vector<atomVector> chainAtomVector;

class Chain {
public:
    Chain(std::vector<std::shared_ptr<Residue>>, std::string);
    int n_residues() { return nResidues; };
    int n_atoms() { return nAtoms; };

    std::vector<std::shared_ptr<Residue>> getResidues() { return residues; }

    chainAtomVector getBackboneAtoms() {
        std::vector<std::vector<std::shared_ptr<Atom>>> backboneAtoms;
        backboneAtoms.reserve((nResidues));
        for(auto const& residue: residues) {
            backboneAtoms.push_back(residue->getBackbone());
        }
        return backboneAtoms;
    }

private:
    std::string chainName;
    std::vector<std::shared_ptr<Residue>> residues;
    int nResidues;
    int nAtoms;
};


#endif //PROSTRUCT_CHAIN_H
