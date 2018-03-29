//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_CHAIN_H
#define PROSTRUCT_CHAIN_H

#include "residue.h"

class Chain {
public:
    Chain(std::vector<std::shared_ptr<Residue>>, std::string);

private:
    std::string chainName;
    int nResidues;
};


#endif //PROSTRUCT_CHAIN_H
