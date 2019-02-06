//
// Created by gil on 27/03/18.
//

#include "prostruct/pdb/PDB.h"

template <typename T>
PDB<T>::PDB(std::string filename_) {

    filename = filename_;

    std::map<std::string, std::map<std::string, atomVector<T>, AASequenceOrder>> chainAtomMap;

    createMap(filename, chainAtomMap, chainOrder);

    nAtoms = 0;
    nResidues = 0;

    xyz.resize(3, 0);

    for (const auto& chain: chainOrder) {

        residueVector<T> residues;

        for (auto const &atomPair: chainAtomMap[chain]) {

            try {
                auto residue = std::make_shared<Residue<T>>(atomPair.second, atomPair.first.substr(0, 3), atomPair.first);
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
            chainMap[chain] = std::make_shared<Chain<T>>(residues, chain);
            nResidues += static_cast<arma::uword>(residues.size());
        }
        catch(const char* msg){
            std::cout << "Chain: " << chain << std::endl;
            std::cout << msg << std::endl;
        }
    }

    numberOfChains = static_cast<int>(chainMap.size());
}

template <typename T>
PDB<T> PDB<T>::fetch(std::string PDB_id) {

    throw "Not implemented";

}

template <typename T>
void PDB<T>::getBackboneAtoms(arma::Mat<T> &C_coords, arma::Mat<T> &O_coords, arma::Mat<T> &N_coords, arma::Mat<T> &CA_coords) {

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

template <typename T>
void PDB<T>::internalKS (arma::Mat<T> &E) {

    arma::Mat<T> C_coords(3, nResidues);
    arma::Mat<T> O_coords(3, nResidues);
    arma::Mat<T> N_coords(3, nResidues);
    arma::Mat<T> CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    std::vector<bool> hasHbond(nResidues, false);

    kabsch_sander(C_coords, O_coords, N_coords, CA_coords, hasHbond, E, nResidues);
}

template <typename T>
arma::Mat<T> PDB<T>::calculate_KabschSander() {

    arma::Mat<T> E(nResidues, nResidues);

    internalKS(E);

    return E;
};

template <typename T>
arma::Col<T> PDB<T>::calculate_ASA(T probe) {

    // calculates atom surface accessibility using the Shrake-Rupley algorithm

    arma::Col<T> asa(static_cast<const arma::uword>(nAtoms));

    shrake_rupley(xyz, radii, asa, nAtoms, probe);

    return asa;
}

template <typename T>
arma::Mat<T> PDB<T>::predict_backboneHbonds() {

    arma::Mat<T> E(nResidues, nResidues);

    internalKS(E);

    E.for_each([](typename arma::Mat<T>::elem_type& elem) {elem = elem < -0.5;});

    return E;
}

template <typename T>
void PDB<T>::calculate_dssp() {

    arma::Mat<T> C_coords(3, nResidues);
    arma::Mat<T> O_coords(3, nResidues);
    arma::Mat<T> N_coords(3, nResidues);
    arma::Mat<T> CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    dssp(C_coords, O_coords, N_coords, CA_coords);
}

template <typename T>
double PDB<T>::calculate_RMSD(PDB &other) {

    // first check if the size is the same
    if (nAtoms != other.n_atoms()) {
        throw "Atom number mismatch";
    }

    return rmsd(xyz, other.getXYZ());
}

template <typename T>
arma::Col<T> PDB<T>::calculate_centroid() {

    arma::Col<T> result(3);

    get_centroid(xyz, result);

    return result;
}

template <typename T>
void PDB<T>::recentre() {

    recentre_molecule(xyz);

}

template <typename T>
arma::Mat<T> PDB<T>::calculate_phi_psi() {

    // calculate phi and psi dihedral angles along the protein backbone
    // returns a matrix of shape (2, n_atoms):
    //  - row 0: phi
    //  - row 1: psi

    // different chains are not connected, therefore need ignore the C and N terminus of consecutive chains

    arma::Mat<T> result(nResidues-chainOrder.size()+1, 2, arma::fill::zeros);

    arma::uword atomPosition = 0;

    arma::Mat<T> C_coords(3, nResidues);
    arma::Mat<T> O_coords(3, nResidues);
    arma::Mat<T> N_coords(3, nResidues);
    arma::Mat<T> CA_coords(3, nResidues);

    getBackboneAtoms(C_coords, O_coords, N_coords, CA_coords);

    for (const auto &chainName: chainOrder) {
        std::cout << atomPosition;
        auto chain = chainMap[chainName];
        arma::uword n_residues = chain->n_residues();

        std::cout << n_residues;

        // create a 3D tensor with all the relevant coordinates
        arma::Cube<T> phiAtomCoords(3,4, n_residues-1);
        arma::Cube<T> psiAtomCoords(3,4, n_residues-1);

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

        arma::Col<T> phi(n_residues-1);
        arma::Col<T> psi(n_residues-1);

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

//void rotate(arma::Col<T> &rotation) {
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

template <typename T>
void PDB<T>::kabsch_rotation(PDB<T> &other) {

    // make copy of xyz
    auto xyz_copy = other.getXYZ();

    kabsch_rotation_(xyz, xyz_copy);

}

template <typename T>
T PDB<T>::kabsch_rmsd(PDB &other) {

    // make copy of xyz
    auto xyz_copy = xyz;
    auto xyz_other_copy = other.getXYZ();

    return kabsch_rmsd_(xyz_copy, xyz_other_copy);

}

template class PDB<float>;
template class PDB<double>;