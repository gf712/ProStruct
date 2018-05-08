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
                radii.insert_rows(nAtoms, residue->getRadii());
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

    arma::vec asa(nAtoms);

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

    arma::mat result(nResidues, 2);

    arma::mat C_coords(3, nResidues);
    arma::mat O_coords(3, nResidues);
    arma::mat N_coords(3, nResidues);
    arma::mat CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    arma::cube phiAtomCoords(3,4,nResidues-1);
    arma::cube psiAtomCoords(3,4,nResidues-1);

    // create a 3D tensor with all the relevant coordinates

    arma::mat &phi_slice_ref = phiAtomCoords.slice(0);
    phi_slice_ref(arma::span::all, 0) = C_coords.col(0);
    phi_slice_ref(arma::span::all, 1) = N_coords.col(0);
    phi_slice_ref(arma::span::all, 2) = CA_coords.col(1);
    phi_slice_ref(arma::span::all, 3) = C_coords.col(1);

    arma::mat &psi_slice_ref = psiAtomCoords.slice(nResidues-1);
    psi_slice_ref(arma::span::all, 0) = N_coords.col(nResidues-1);
    psi_slice_ref(arma::span::all, 1) = CA_coords.col(nResidues);
    psi_slice_ref(arma::span::all, 2) = C_coords.col(nResidues);
    psi_slice_ref(arma::span::all, 3) = N_coords.col(nResidues);

    for (int i = 1; i < nResidues-1; ++i) {
        // phi angle
        phi_slice_ref = phiAtomCoords.slice(i);
        phi_slice_ref(arma::span::all, 0) = C_coords.col(i);
        phi_slice_ref(arma::span::all, 1) = N_coords.col(i);
        phi_slice_ref(arma::span::all, 2) = CA_coords.col(i+1);
        phi_slice_ref(arma::span::all, 3) = C_coords.col(i+1);

        // psi angle
        psi_slice_ref = psiAtomCoords.slice(i);
        psi_slice_ref(arma::span::all, 0) = N_coords.col(i);
        psi_slice_ref(arma::span::all, 1) = CA_coords.col(i+1);
        psi_slice_ref(arma::span::all, 2) = C_coords.col(i+1);
        psi_slice_ref(arma::span::all, 3) = N_coords.col(i+1);

    }

    arma::vec phi(nResidues-1);
    arma::vec psi(nResidues-1);

    dihedrals(phiAtomCoords, phi);
    dihedrals(psiAtomCoords, psi);

    result(arma::span(0, nResidues-1), 0) = phi;
    result(arma::span(1, nResidues), 0) = psi;
}