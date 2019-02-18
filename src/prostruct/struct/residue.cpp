//
// Created by gil on 27/03/18.
//

#include "prostruct/struct/residue.h"
#include <algorithm>
#include <iostream>

template <typename T>
Residue<T>::Residue(atomVector<T> atoms_, const std::string& aminoAcidName_, const std::string& residueName_)
{

	/**
     *  The Residue class represent one of the twenty standard amino acids.
     *  A Residue is made of Atom objects which are connected via Bond objects.
     *
     *  @param[in] atoms_ Vector of Atom objects
     *  @param[in] aminoAcidName_ name of amino acid
     *  @param[in] residueName_ name of residue
     */

	backbone = std::vector<int>(4);

	if (aminoAcidIndex.find(aminoAcidName_) != aminoAcidIndex.end()) {
		aminoAcidName = aminoAcidName_;
		aminoAcid = static_cast<AminoAcid>(aminoAcidIndex.at(aminoAcidName));
	} else
		throw "Unknown amino acid!";

	// each atom is responsible to form a bond with the previous atom
	int i = 0;
	for (auto atom : atoms_) {

		// assumes that we are using the following scheme:
		// amide group:
		// N -> CA <- C <- O
		// CA -> is bound to the sidechain
		// CA ->Glycine CB -> CG
		// and atoms are passed in the same order as in PDBs:
		// N, CA, C, O, R (CB,..)
		// plus the special rules for the cyclic amino acids

		auto name = atom->getName();
		// is the atom in the backbone?
		auto aaBackbone = backboneIndexMap.find(name);
		auto aaLocation_ = (aaBackbone != backboneIndexMap.end()) ? aaLocation::Backbone : aaLocation::Sidechain;

		switch (aaLocation_) {
		case aaLocation::Backbone: {
			// inserts backbone atom in correct location -> this is important because we will always assume
			// that the N is at position 0 of atoms and C at position 2.
			backbone.at(backboneIndexMap.at(name)) = i;
		} break;
		case aaLocation::Sidechain:
			sidechain.push_back(i);
		}

		atomMap[name] = i;
		atoms.emplace_back(atom);
		atom->setRadius(aminoAcidRadii.at(static_cast<int>(aminoAcid)).at(name));
		i++;
	}

	if (backbone.size() != 4) {
		throw "Expected four atoms in the backbone, got " + std::to_string(backbone.size());
	}

	residueName = residueName_;

	createBonds();
}

template <typename T>
void Residue<T>::createBonds()
{

	bool first = true;
	std::vector<std::string> atomPair;
	auto res = aminoAcidAtoms[static_cast<int>(aminoAcid)];

	std::vector<int> positions;
	std::copy(backbone.begin(), backbone.end(), std::back_inserter(positions));
	std::copy(sidechain.begin(), sidechain.end(), std::back_inserter(positions));

	for (auto const& pos : positions) {

		auto atom = atoms[pos];

		auto name = atom->getName();

		if (first) {
			first = false;
		}

		else {

			// checks if the atom is expected
			if (res.find(name) != res.end()) {
				atomPair = { name, res[name] };
			} else
				throw "Unknown atom: " + name;

			//            std::cout << "Bond: from: " << atomPair[0] << " to: " << atomPair[1] << std::endl;

			// does bond exist?
			if (!atom->hasBond(atoms[atomMap[atomPair[1]]])) {
				atom->addBond(atoms[atomMap[atomPair[1]]], 1);
				//                std::cout << atom->getName() << " n bond: " << atom->getNumberOfBonds()
				//                          << " " << atoms[atomMap[atomPair[1]]]->getName() << " n bonds: "
				//                          << atoms[atomMap[atomPair[1]]]->getNumberOfBonds() << std::endl;
			}
		}
		//        std::cout << "Completed " << atom->getName() << std::endl;
	}

	switch (aminoAcid) {

	case AminoAcid::PRO:
		atoms[atomMap["CD"]]->addBond(atoms[atomMap["N"]], 1);
		break;
	case AminoAcid::TRP: {
		atoms[atomMap["CD2"]]->addBond(atoms[atomMap["CE2"]], 1);
		atoms[atomMap["CH2"]]->addBond(atoms[atomMap["CZ3"]], 1);
	} break;
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

	xyz.resize(3, atoms.size());
	radii.resize(atoms.size());

	arma::uword i = 0;
	for (const auto& atom : getBackbone()) {
		xyz.at(0, i) = atom->getX();
		xyz.at(1, i) = atom->getY();
		xyz.at(2, i) = atom->getZ();
		radii.at(i) = atom->getRadius();
		i++;
	}
	for (const auto& atom : getSidechain()) {
		xyz.at(0, i) = atom->getX();
		xyz.at(1, i) = atom->getY();
		xyz.at(2, i) = atom->getZ();
		radii.at(i) = atom->getRadius();
		i++;
	}
}

template <typename T>
void Residue<T>::link(std::shared_ptr<Residue<T>> residue_)
{
	// links this (C-terminus) with residue_ (N-terminus)
	atoms[backbone[0]]->addBond(residue_->getBackbone()[2], 1);
}

template class Residue<float>;
template class Residue<double>;
