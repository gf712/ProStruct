//
// Created by gil on 27/03/18.
//

#include "../include/PDB.h"

PDB::PDB(std::string filename_) {

    filename = filename_;

    std::map<std::string, std::map<std::string, std::vector<std::shared_ptr<Atom>>, AASequenceOrder>> chainAtomMap;

    createMap(filename, chainAtomMap, chainOrder);

    nAtoms = 0;
    nResidues = 0;

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
            nResidues += static_cast<arma::uword>(residues.size());
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

void PDB::getBackboneAtoms(arma::mat &C_coords, arma::mat &O_coords, arma::mat &N_coords, arma::mat &CA_coords) {

    // get all C, O, N atoms positions
    arma::uword i = 0;
    arma::uword pos = 0;

    for (auto const& chain: chainOrder) {

        for (auto const& residue: chainMap[chain]->getResidues()) {

            for (arma::uword coord = 0; coord < 3; ++coord) {

                N_coords.at(coord, i) = xyz.at(coord, pos);
                CA_coords.at(coord, i) = xyz.at(coord, pos + 1);
                C_coords.at(coord, i) = xyz.at(coord, pos + 2);
                O_coords.at(coord, i) = xyz.at(coord, pos + 3);

            }

            pos+=residue->n_atoms();
            i++;
        }
    }
}


void PDB::internalKS (arma::mat &E) {

    arma::mat C_coords(3, nResidues);
    arma::mat O_coords(3, nResidues);
    arma::mat N_coords(3, nResidues);
    arma::mat CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    std::vector<bool> hasHbond(nResidues, false);

    kabsch_sander(C_coords, O_coords, N_coords, CA_coords, hasHbond, E, nResidues);
}


arma::mat PDB::calculate_KabschSander() {

    arma::mat E(nResidues, nResidues);

    internalKS(E);

    return E;
};


arma::mat PDB::predict_backboneHbonds() {

    arma::mat E(nResidues, nResidues);

    internalKS(E);

    E.for_each([] (arma::mat::elem_type& elem) {elem = elem < -0.5;});

    return E;
}


void PDB::calculate_dssp() {

    int n_residues = 0;

    for (auto const& chain: chainMap) {

        n_residues += static_cast<arma::uword>(chain.second->n_residues());

    }

    arma::mat C_coords(3, n_residues);
    arma::mat O_coords(3, n_residues);
    arma::mat N_coords(3, n_residues);
    arma::mat CA_coords(3, n_residues);

    // get all C, O, N atoms positions
    arma::uword i = 0;
    arma::uword pos = 0;

    for (auto const& chain: chainOrder) {

        for (auto const& residue: chainMap[chain]->getResidues()) {

            for (arma::uword coord = 0; coord < 3; ++coord) {

                N_coords.at(coord, i) = xyz.at(coord, pos);
                CA_coords.at(coord, i) = xyz.at(coord, pos + 1);
                C_coords.at(coord, i) = xyz.at(coord, pos + 2);
                O_coords.at(coord, i) = xyz.at(coord, pos + 3);

            }

            pos+=residue->n_atoms();
            i++;
        }
    }

    dssp(C_coords, O_coords, N_coords, CA_coords);
}
