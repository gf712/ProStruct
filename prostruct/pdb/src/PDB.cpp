//
// Created by gil on 27/03/18.
//

#include "../include/PDB.h"
#include "PDBparser.h"
#include <iostream>

PDB::PDB(std::string filename_) {

    filename = filename_;
    std::map<std::string, std::map<std::string, std::vector<std::shared_ptr<Atom>>, AASequenceOrder>> ChainMap;

    createMap(filename, ChainMap);

    numberOfChains = static_cast<int>(ChainMap.size());

    for (const auto& pair: ChainMap) {

        std::vector<std::shared_ptr<Residue>> residues;


        for (auto const &atomPair: pair.second) {

            try {
                residues.emplace_back(std::make_shared<Residue>(atomPair.second, atomPair.first.substr(0, 3), atomPair.first));
            }
            catch(const char* msg){
                std::cout << "Residue: " << atomPair.first << std::endl;
                std::cout << msg << std::endl;
            }

        }

        try {
            auto chain = std::make_shared<Chain>(residues, pair.first);
        }
        catch(const char* msg){
            std::cout << "Chain: " << pair.first << std::endl;
            std::cout << msg << std::endl;
        }
    }
}

PDB PDB::fetch(std::string PDB_id) {

    throw "Not impletemented";

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