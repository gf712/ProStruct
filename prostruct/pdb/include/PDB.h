//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#include <string>
#include <stdio.h>
#include "chain.h"

class PDB {

public:

    PDB(std::string filename);

    ~PDB() = default;

    static PDB fetch(std::string);

    std::string getFilename() { return filename;}

    int n_chains() { return numberOfChains; }

    std::vector<std::string> getChainIDs() {
        std::vector<std::string> chainIDs;
        for (auto const& pair: chainMap) {
            chainIDs.push_back(pair.first);
        }
        return chainIDs;
    }

    std::shared_ptr<Chain> getChain(std::string name_) { return chainMap[name_]; }

private:

    std::string filename;
    int numberOfChains;
    std::map<std::string, std::shared_ptr<Chain>> chainMap;

};


#endif //PROSTRUCT_PDB_H
