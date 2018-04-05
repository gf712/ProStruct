//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include "residue.h"

class Chain {
public:
    Chain(std::vector<std::shared_ptr<Residue>>, std::string);
    int n_residues() { return nResidues; };
    int n_atoms() { return nAtoms; };

private:
    std::string chainName;
    std::vector<std::shared_ptr<Residue>> residues;
    int nResidues;
    int nAtoms;
};


#endif //PROSTRUCT_CHAIN_H
