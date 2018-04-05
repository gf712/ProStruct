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

    // Overload operator
//    friend std::istream& operator>>(std::istream& str, std::shared_ptr<Chain>& chain)
//    {
//        return nullptr;
//    }


private:

//    void parse(FILE);

    std::string filename;
    std::string title;

};


#endif //PROSTRUCT_PDB_H
