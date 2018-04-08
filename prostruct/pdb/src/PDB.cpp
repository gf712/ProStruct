//
// Created by gil on 27/03/18.
//

#include "../include/PDB.h"

PDB::PDB(std::string filename_) {

    filename = filename_;

    std::map<std::string, std::map<std::string, std::vector<std::shared_ptr<Atom>>, AASequenceOrder>> chainAtomMap;

    createMap(filename, chainAtomMap, chainOrder);

    nAtoms = 0;

    xyz.resize(3, 0);

    for (const auto& chain: chainOrder) {

        std::vector<std::shared_ptr<Residue>> residues;

        for (auto const &atomPair: chainAtomMap[chain]) {

            try {
                auto residue = std::make_shared<Residue>(atomPair.second, atomPair.first.substr(0, 3), atomPair.first);
                residues.emplace_back(residue);
                xyz.insert_cols(nAtoms, residue->getXYZ());
                nAtoms += residue->n_atoms();
//                residue->getXYZ().print();
            }
            catch(const char* msg){
                std::cout << "Residue: " << atomPair.first << std::endl;
                std::cout << msg << std::endl;
            }
        }

        try {
            chainMap[chain] = std::make_shared<Chain>(residues, chain);
        }
        catch(const char* msg){
            std::cout << "Chain: " << chain << std::endl;
            std::cout << msg << std::endl;
        }
    }

    numberOfChains = static_cast<int>(chainMap.size());
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