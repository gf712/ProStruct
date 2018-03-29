//
// Created by gil on 27/03/18.
//

#include <algorithm>
#include <iostream>
#include "../include/residue.h"


std::map<std::string, int> backboneMap
        {
                {"N",  0},
                {"CA", 1},
                {"C",  2},
                {"O",  3},
        };

// hardcoded amino acid atoms
std::vector<std::map<std::string, std::string>> aminoAcidAtoms
        {
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD", "CG"   }, {"CD", "HD2"  }, {"CD", "HD3"},   {"NE", "CD"},{"CG", "HG2"},{"CG", "HG3"},{"CZ", "NE"},{"NH1", "CZ"},{"NH2", "CZ"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HE", "NE"},{"HH11", "NH1"},{"HH12", "NH1"},{"HH21", "NH2"},{"HH22", "NH2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"HB1","CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"H", "N"     }, {"H2", "N"    }, {"H3", "N"},     {"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"ND2", "CG"  }, {"OD1", "CG"  }, {"H", "N"},      {"H2", "N"},{"H3", "N"},{"HD21", "ND2"},{"HD22", "ND2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"OD1", "CG"  }, {"OD2", "CG"  }, {"H", "N"},      {"H2", "N"},{"H3", "N"},{"HD2", "OD2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"HB2","CB"}, {"HB3", "CB"}, {"SG" , "CB"}, {"H", "N"     }, {"H2", "N"    }, {"H3", "N"},     {"HG", "SG"},{"HXT", "OXT"},},
                {{"-C", "N"}, {"C", "CA"}, {"O", "C"}, {"C", "OXT"}, {"CB", "CA"}, {"CA", "HA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD", "CG"   }, {"NE2", "CD"  }, {"OE1", "CD"},{"CG", "HG2"},{"CG", "HG3"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HE21", "NE2"},{"HE22", "NE2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD", "CG"   }, {"OE1", "CD"  }, {"OE2", "CD"},   {"CG", "HG2"},{"CG", "HG3"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HE2", "OE2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"HA2","CA"}, {"HA3","CA"}, {"CA", "N"}, {"H",  "N" }, {"H2",  "N" }, {"H3" , "N" }, {"HXT", "OXT" }},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD2", "CG"  }, {"HD2", "CD2" }, {"HE1", "CE1"},{"CE1", "ND1"},{"NE2", "CE1"},{"ND1", "CG"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HD1", "ND1"},{"HE2", "NE2"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG1","CB"}, {"CG2", "CB"}, {"HB" , "CB"}, {"CD1", "CG1" }, {"CD1", "HD11"}, {"CD1", "HD12"}, {"CD1", "HD13"},{"CG1", "HG12"},{"CG1", "HG13"},{"CG2", "HG21"},{"CG2", "HG22"},{"CG2", "HG23"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD1", "CG"  }, {"CD1", "HD11"}, {"CD1", "HD12"}, {"CD1", "HD13"},{"CD2", "CG"},{"CD2", "HD21"},{"CD2", "HD22"},{"CD2", "HD23"},{"HG","CG"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CE", "CD"   }, {"CD", "CG"   }, {"CD", "HD2"},   {"CD", "HD3"},{"CE", "HE2"},{"CE", "HE3"},{"NZ", "CE"},{"HG2", "CG"},{"HG3","CG"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"},{"HZ1", "NZ"},{"HZ2", "NZ"},{"HZ3", "NZ"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"HE1", "CE"  }, {"HE2", "CE"  }, {"HE3", "CE"},   {"SD", "CG"}, {"CE", "SD"},{"CG", "HG2"},{"CG", "HG3"},{"CG", "SD"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CE1", "CD1" }, {"CD1", "CG"  }, {"HD1", "CD1"},  {"CE2", "CD2"},{"CD2", "CG"},{"HD2", "CD2"},{"CE1", "CZ"},{"CZ", "CE2"},{"HE2", "CE2"},{"HZ", "CZ"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD", "CG"   }, {"CD", "HD2"  }, {"CD", "HD3"},   {"CD", "N"},{"CG", "HG2"},{"CG", "HG3"},{"H", "N"},{"H3", "N"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"HB2","CB"}, {"HB3", "CB"}, {"OG" , "CB"}, {"H", "N"     }, {"H2", "N"    }, {"H3", "N"},     {"HG", "OG"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG2","CB"}, {"HB" , "CB"}, {"OG1", "CB"}, {"HG21", "CG2"}, {"HG22", "CG2"}, {"HG23", "CG2"}, {"H", "N"},{"H2", "N"},{"H3", "N"},{"HG1", "OG1"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CD1", "CG"  }, {"CD1", "HD1" }, {"NE1", "CD1"},  {"CE3", "CD2"},{"CD2", "CG"},{"CZ2", "CE2"},{"CE2", "NE1"},{"CZ3", "CE3"},{"HE3", "CE3"},{"CH2", "CZ2"},{"CH2", "CZ3"},{"HH2", "CH2"},{"HZ2", "CZ2"},{"HZ3", "CZ3"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HE1", "NE1"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG", "CB"}, {"HB2", "CB"}, {"HB3", "CB"}, {"CE1", "CD1" }, {"CD1", "CG"  }, {"HD1", "CD1"},  {"CE2", "CD2"}, {"CD2", "CG"},  {"CD2", "HD2"}, {"CE1", "CD1"}, {"CZ", "CE1"},{"HE1", "CE1"},{"CE2", "CD2"},{"HE2", "CE2"},{"OH", "CZ"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HH", "OH"},{"HXT", "OXT"},},
                {{"N", "-C"}, {"C", "CA"}, {"O", "C"}, {"OXT", "C"}, {"CB", "CA"}, {"HA", "CA"}, {"CA", "N"}, {"CG1","CB"}, {"CG2", "CB"}, {"HB" , "CB"}, {"CG1", "HG11"}, {"CG1", "HG12"}, {"CG1", "HG13"}, {"CG2", "HG21"},{"CG2", "HG22"},{"CG2", "HG23"},{"H", "N"},{"H2", "N"},{"H3", "N"},{"HXT", "OXT"}}
        };

// hardcoded map with amino acid atoms
std::map<std::string, int> aminoAcidIndex
        {
                {"ARG",  0},
                {"ALA",  1},
                {"ASN",  2},
                {"ASP",  3},
                {"CYS",  4},
                {"GLN",  5},
                {"GLU",  6},
                {"GLY",  7},
                {"HIS",  8},
                {"ILE",  9},
                {"LEU", 10},
                {"LYS", 11},
                {"MET", 12},
                {"PHE", 13},
                {"PRO", 14},
                {"SER", 15},
                {"THR", 16},
                {"TRP", 17},
                {"TYR", 18},
                {"VAL", 19}
        };


Residue::Residue(std::vector<std::shared_ptr<Atom>> atoms_, std::string aminoAcidName_, std::string residueName_) {

    /**
     *  The Residue class represent one of the twenty standard amino acids.
     *  A Residue is made of Atom objects which are connected via Bond objects.
     *
     *  @param[in] atoms_ Vector of Atom objects
     *  @param[in] aminoAcidName_ name of amino acid
     *  @param[in] residueName_ name of residue
     */

    backbone = std::vector<int>(4);

    // each atom is responsible to form a bond with the previous atom
    int i = 0;
    for (auto atom: atoms_) {

        // assumes that we are using the following scheme:
        // amide group:
        // N -> CA <- C <- O
        // CA -> is bound to the sidechain
        // CA -> CB -> CG
        // and atoms are passed in the same order as in PDBs:
        // N, CA, C, O, R (CB,..)
        // plus the special rules for the cyclic amino acids

        auto name = atom->getName();
        // is the atom in the backbone?
        auto aaBackbone = backboneMap.find(name);
        auto aaLocation_ = (aaBackbone != backboneMap.end()) ? aaLocation::Backbone : aaLocation::Sidechain;

        switch (aaLocation_) {
            case aaLocation::Backbone: {
                // inserts backbone atom in correct location -> this is important because we will always assume
                // that the N is at position 0 of atoms and C at position 2.
                backbone[backboneMap[name]] = i;
            }
                break;
            case aaLocation::Sidechain:
                sidechain.push_back(i);
        }

        atomMap[atom->getName()] = i;
        atoms.emplace_back(atom);
        i++;
    }

    if (backbone.size() != 4) {
        throw "Expected four atoms in the backbone, got " + std::to_string(backbone.size());
    }

    if (aminoAcidIndex.find(aminoAcidName_) != aminoAcidIndex.end()) {
        aminoAcidName = aminoAcidName_;
        aminoAcid = static_cast<AminoAcid>(aminoAcidIndex[aminoAcidName]);
    }
    else
        throw "Unknown amino acid!";

    residueName = residueName_;

    createBonds();

}

void Residue::createBonds() {

    bool first = true;
    std::vector<std::string> atomPair;
    auto res = aminoAcidAtoms[static_cast<int>(aminoAcid)];

    std::vector<int> positions;
    std::copy(backbone.begin(), backbone.end(), std::back_inserter(positions));
    std::copy(sidechain.begin(), sidechain.end(), std::back_inserter(positions));

    for (auto const &pos: positions) {

        auto atom = atoms[pos];

        auto name = atom->getName();

        std::cout << "Name: " << atom->getName() << std::endl;

        if (first) {
            first=false;
        }

        else {

            // checks if the atom is expected
            if (res.find(name) != res.end()) {
                atomPair = {name, res[name]};
            }
            else throw "Unknown atom!";

            std::cout << "Bond: from: " << atomPair[0] << " to: " << atomPair[1] << std::endl;

            // does bond exist?
            if (!atom->hasBond(atoms[atomMap[atomPair[1]]])) {
                atom->addBond(atoms[atomMap[atomPair[1]]], 1);
                std::cout << atom->getName() << " n bond: " << atom->getNumberOfBonds()
                          << " " << atoms[atomMap[atomPair[1]]]->getName() << " n bonds: "
                          << atoms[atomMap[atomPair[1]]]->getNumberOfBonds() << std::endl;
            }
        }
        std::cout << "Completed " << atom->getName() << std::endl;
    }

    switch (aminoAcid) {

        case AminoAcid::PRO:
            atoms[atomMap["CD"]]->addBond(atoms[atomMap["N"]], 1);
            break;
        case AminoAcid::TRP: {
            atoms[atomMap["CD2"]]->addBond(atoms[atomMap["CE2"]], 1);
            atoms[atomMap["CH2"]]->addBond(atoms[atomMap["CZ3"]], 1);
        }
            break;
        case AminoAcid::HIS:
            atoms[atomMap["CD2"]]->addBond(atoms[atomMap["NE2"]], 1);
            break;
        case AminoAcid::PHE:
            atoms[atomMap["CE1"]]->addBond(atoms[atomMap["CZ"]], 1);
            break;
        case AminoAcid::TYR:
            atoms[atomMap["CE2"]]->addBond(atoms[atomMap["OH"]], 1);
        default:
            break;
    }

    std::cout << "Completed bond creation!" << std::endl;
}

void Residue::link(std::shared_ptr<Residue> residue_) {
    // links this (C-terminus) with residue_ (N-terminus)
    atoms[backbone[0]]->addBond(residue_->getBackbone()[2], 1);
}
