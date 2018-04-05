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

private:

    std::string filename;
    std::string title;
    int numberOfChains;

};


#endif //PROSTRUCT_PDB_H
