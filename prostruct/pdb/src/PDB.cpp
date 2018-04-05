//
// Created by gil on 27/03/18.
//

#include "../include/PDB.h"
#include "PDBparser.h"


PDB::PDB(std::string filename_) {

    filename = filename_;
    std::map<std::string, std::map<std::string, std::vector<std::shared_ptr<Atom>>, AASequenceOrder>> ChainMap;

    createMap(filename, ChainMap);

}

PDB PDB::fetch(std::string PDB_id) {



}

//
//void PDB::parse(FILE file) {
//
//    std::string line;
//
//    while(std::getline(str, file)) {
//
//    }
//
//}