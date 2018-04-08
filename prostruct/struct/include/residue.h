//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_RESIDUE_H
#define PROSTRUCT_RESIDUE_H

#include "atom.h"
#include <armadillo>

enum class aaLocation {
    Backbone,
    Sidechain
};

enum class AminoAcid {
    ARG,
    ALA,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    ILE,
    LEU,
    LYS,
    MET,
    PHE,
    PRO,
    SER,
    THR,
    TRP,
    TYR,
    VAL
};

class Residue {
public:
    // takes an arbitrary number of atoms and tries to form a residue
//    Residue(std::unique_ptr<Atom> atoms...);
    Residue(std::vector<std::shared_ptr<Atom>>, std::string, std::string);

    std::vector<std::shared_ptr<Atom>> getBackbone() {
        std::vector<std::shared_ptr<Atom>> backboneAtoms;
        for (auto const &i: backbone) {
            backboneAtoms.emplace_back(atoms[i]);
        }
        return backboneAtoms;
    }
    std::vector<std::shared_ptr<Atom>> getSidechain() {
        std::vector<std::shared_ptr<Atom>> sidechainAtoms;
        for (auto const &i: sidechain) {
            sidechainAtoms.emplace_back(atoms[i]);
        }
        return sidechainAtoms;
    }

    std::shared_ptr<Atom> operator[](const int index) { return atoms[index]; }

    arma::mat getXYZ() { return xyz; }

    std::shared_ptr<Atom> getAtom(int index) { return atoms[0]; }

    std::string getResidueName() { return residueName; }

    void link(std::shared_ptr<Residue>);

    void createBonds();

    int n_atoms() { return atoms.size(); };

    void print_atoms() {

        for (auto const &atom: atoms) {
            std::cout << atom->getName() << "-";
        }
    }

private:
    arma::mat xyz;
    std::vector<int> backbone; /**< A vector with the index number of the backbone atoms */
    std::vector<int> sidechain; /**< A vector with the index number of the sidechain atoms */
    std::string aminoAcidName; /**< Name of the amino acid, e.g. ALA */
    enum AminoAcid aminoAcid; /**< Amino acid enum, e.g. ALA */
    std::string residueName; /**< Name of the residue, e.g. ALA1 */
    std::vector<std::shared_ptr<Atom>> atoms; /**< A vector with the pointers to the Atom objects */
    std::map<std::string, int> atomMap; /**< Map atom name to internal index */

};


#endif //PROSTRUCT_RESIDUE_H
