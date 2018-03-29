//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#include <string>

class PDB {

public:
    PDB(std::string filename) {
        filename = filename;
    }
    ~PDB() = default;


private:
    std::string filename;
    std::string title;
};


#endif //PROSTRUCT_PDB_H
