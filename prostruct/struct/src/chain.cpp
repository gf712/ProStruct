//
// Created by gil on 27/03/18.
//

#include "../include/chain.h"

Chain::Chain(std::vector<std::shared_ptr<Residue>> residues_, std::string chainName_) {

    bool first = true;
    std::shared_ptr<Residue> previousResidue;

    // link residues
    for (auto const &residue: residues_) {

        if (first) {
            previousResidue = residue;
            first = false;
        }

        else {
            residue->link(previousResidue);
        }
    }

    chainName = chainName_;
    nResidues = residues_.size();

}