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
                xyz.insert_cols(static_cast<const arma::uword>(nAtoms), residue->getXYZ());
                radii.insert_rows(static_cast<const arma::uword>(nAtoms), residue->getRadii());
                nAtoms += residue->n_atoms();
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

    throw "Not implemented";

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


arma::vec PDB::calculate_ASA(double probe) {

    // calculates atom surface accessibility using the Shrake-Rupley algorithm

    arma::vec asa(static_cast<const arma::uword>(nAtoms));

    shrake_rupley(xyz, radii, asa, nAtoms, probe);

    return asa;
}


arma::mat PDB::predict_backboneHbonds() {

    arma::mat E(nResidues, nResidues);

    internalKS(E);

    E.for_each([] (arma::mat::elem_type& elem) {elem = elem < -0.5;});

    return E;
}


void PDB::calculate_dssp() {

    arma::mat C_coords(3, nResidues);
    arma::mat O_coords(3, nResidues);
    arma::mat N_coords(3, nResidues);
    arma::mat CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    dssp(C_coords, O_coords, N_coords, CA_coords);
}

double PDB::calculate_RMSD(PDB &other) {

    // first check if the size is the same
    if (nAtoms != other.n_atoms()) {
        throw "Atom number mismatch";
    }

    return rmsd(xyz, other.getXYZ());
}


arma::vec PDB::calculate_centroid() {

    arma::vec result(3);

    get_centroid(xyz, result);

    return result;
}

void PDB::recentre() {

    recentre_molecule(xyz);

}

arma::mat PDB::calculate_phi_psi() {

    // calculate phi and psi dihedral angles along the protein backbone
    // returns a matrix of shape (2, n_atoms):
    //  - row 0: phi
    //  - row 1: psi

    // different chains are not connected, therefore need ignore the C and N terminus of consecutive chains

    arma::mat result(nResidues-chainOrder.size()+1, 2, arma::fill::zeros);

    arma::uword atomPosition = 0;

    arma::mat C_coords(3, nResidues);
    arma::mat O_coords(3, nResidues);
    arma::mat N_coords(3, nResidues);
    arma::mat CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    for (const auto &chainName: chainOrder) {
        std::cout << atomPosition;
        auto chain = chainMap[chainName];
        auto n_residues = static_cast<arma::uword>(chain->n_residues());

        std::cout << n_residues;

        // create a 3D tensor with all the relevant coordinates
        arma::cube phiAtomCoords(3,4, n_residues-1);
        arma::cube psiAtomCoords(3,4, n_residues-1);

        for (arma::uword i=atomPosition; i<n_residues-1; ++i) {

            // phi angle
            phiAtomCoords.slice(i)(arma::span::all, 0) = C_coords.col(i);
            phiAtomCoords.slice(i)(arma::span::all, 1) = N_coords.col(i+1);
            phiAtomCoords.slice(i)(arma::span::all, 2) = CA_coords.col(i+1);
            phiAtomCoords.slice(i)(arma::span::all, 3) = C_coords.col(i+1);

            // psi angle
            psiAtomCoords.slice(i)(arma::span::all, 0) = N_coords.col(i);
            psiAtomCoords.slice(i)(arma::span::all, 1) = CA_coords.col(i);
            psiAtomCoords.slice(i)(arma::span::all, 2) = C_coords.col(i);
            psiAtomCoords.slice(i)(arma::span::all, 3) = N_coords.col(i+1);

        }

        arma::vec phi(n_residues-1);
        arma::vec psi(n_residues-1);

        std::cout << "Calculating phi/psi" << std::endl;

        dihedrals(phiAtomCoords, phi);
        dihedrals(psiAtomCoords, psi);

        (phi * (180.0 / M_PI)).print();

        result(arma::span(atomPosition, n_residues-2), 0) = phi * (180.0 / M_PI);
        result(arma::span(atomPosition+1, n_residues-1), 1) = psi * (180.0 / M_PI);

        atomPosition+=n_residues+1;

        }

    return result;
}

//void rotate(arma::vec &rotation) {
//
//    // first check which axes need to be rotated
//    arma::uvec mask(3);
//    for (int i = 0; i < 3; ++i) {
//        mask[i] = rotation[i] != 0 ? 1 : 0;
//    }
//
//    // then calculate the rotation matrices
//    for (int i = 0; i < 3; ++i) {
//
//    }
//
//}

void PDB::kabsch_rotation(PDB &other) {

    // make copy of xyz
    auto xyz_copy = other.getXYZ();

    kabsch_rotation_(xyz, xyz_copy);

}

double PDB::kabsch_rmsd(PDB &other) {

    // make copy of xyz
    auto xyz_copy = xyz;
    auto xyz_other_copy = other.getXYZ();

    return kabsch_rmsd_(xyz_copy, xyz_other_copy);

}