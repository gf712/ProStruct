//
// Created by gil on 27/03/18.
//

#include "prostruct/struct/residue.h"
#include <algorithm>
#include <iostream>

using namespace prostruct;

const static stringIndexMap backboneIndexMap = {
	{ "N", 0 },
	{ "CA", 1 },
	{ "C", 2 },
	{ "O", 3 },
};

const static std::vector<std::map<std::string, double>> aminoAcidRadii {
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD", 1.87 }, { "NE", 1.65 }, { "CZ", 1.76 }, { "NH1", 1.65 }, { "NH2", 1.65 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "OD1", 1.40 }, { "ND2", 1.65 }, { "AD1", 1.65 }, { "AD2", 1.65 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "OD1", 1.40 }, { "OD2", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "SG", 1.85 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD", 1.76 }, { "OE1", 1.40 }, { "NE2", 1.65 }, { "AE1", 1.65 }, { "AE2", 1.65 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD", 1.76 }, { "OE1", 1.40 }, { "OE2", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "ND1", 1.65 }, { "CE1", 1.76 }, { "NE2", 1.65 }, { "CD2", 1.76 }, { "AD1", 1.76 }, { "AE1", 1.76 }, { "AE2", 1.76 }, { "AD2", 1.76 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG2", 1.87 }, { "CG1", 1.87 }, { "CD1", 1.87 }, { "CD", 1.87 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD1", 1.87 }, { "CD2", 1.87 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD", 1.87 }, { "CE", 1.87 }, { "NZ", 1.50 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "SD", 1.85 }, { "CE", 1.87 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "CD1", 1.76 }, { "CE1", 1.76 }, { "CZ", 1.76 }, { "CE2", 1.76 }, { "CD2", 1.76 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.87 }, { "CD", 1.87 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "OG", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG2", 1.87 }, { "OG1", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "CD1", 1.76 }, { "NE1", 1.65 }, { "CE2", 1.76 }, { "CZ2", 1.76 }, { "CH2", 1.76 }, { "CZ3", 1.76 }, { "CE3", 1.76 }, { "CD2", 1.76 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG", 1.76 }, { "CD1", 1.76 }, { "CE1", 1.76 }, { "CZ", 1.76 }, { "CE2", 1.76 }, { "CD2", 1.76 }, { "OH", 1.40 }, { "OXT", 1.4 } },
	{ { "N", 1.65 }, { "CA", 1.87 }, { "C", 1.76 }, { "O", 1.40 }, { "CB", 1.87 }, { "CG1", 1.87 }, { "CG2", 1.87 }, { "OXT", 1.4 } }

};

// hardcoded amino acid atoms
const static aminoAcidAtomMap aminoAcidAtoms {
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD", "CG" },
		{ "CD", "HD2" },
		{ "CD", "HD3" },
		{ "NE", "CD" },
		{ "CG", "HG2" },
		{ "CG", "HG3" },
		{ "CZ", "NE" },
		{ "NH1", "CZ" },
		{ "NH2", "CZ" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HE", "NE" },
		{ "HH11", "NH1" },
		{ "HH12", "NH1" },
		{ "HH21", "NH2" },
		{ "HH22", "NH2" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "HB1", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "ND2", "CG" },
		{ "OD1", "CG" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HD21", "ND2" },
		{ "HD22", "ND2" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "OD1", "CG" },
		{ "OD2", "CG" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HD2", "OD2" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "SG", "CB" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HG", "SG" },
		{ "HXT", "OXT" },
	},
	{
		{ "-C", "N" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "C", "OXT" },
		{ "CB", "CA" },
		{ "CA", "HA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD", "CG" },
		{ "NE2", "CD" },
		{ "OE1", "CD" },
		{ "CG", "HG2" },
		{ "CG", "HG3" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HE21", "NE2" },
		{ "HE22", "NE2" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD", "CG" },
		{ "OE1", "CD" },
		{ "OE2", "CD" },
		{ "CG", "HG2" },
		{ "CG", "HG3" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HE2", "OE2" },
		{ "HXT", "OXT" },
	},
	{ { "N", "-C" }, { "C", "CA" }, { "O", "C" }, { "OXT", "C" }, { "HA2", "CA" }, { "HA3", "CA" }, { "CA", "N" }, { "H", "N" }, { "H2", "N" }, { "H3", "N" }, { "HXT", "OXT" } },
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD2", "CG" },
		{ "HD2", "CD2" },
		{ "HE1", "CE1" },
		{ "CE1", "ND1" },
		{ "NE2", "CE1" },
		{ "ND1", "CG" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HD1", "ND1" },
		{ "HE2", "NE2" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG1", "CB" },
		{ "CG2", "CB" },
		{ "HB", "CB" },
		{ "CD1", "CG1" },
		{ "CD1", "HD11" },
		{ "CD1", "HD12" },
		{ "CD1", "HD13" },
		{ "CG1", "HG12" },
		{ "CG1", "HG13" },
		{ "CG2", "HG21" },
		{ "CG2", "HG22" },
		{ "CG2", "HG23" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD1", "CG" },
		{ "CD1", "HD11" },
		{ "CD1", "HD12" },
		{ "CD1", "HD13" },
		{ "CD2", "CG" },
		{ "CD2", "HD21" },
		{ "CD2", "HD22" },
		{ "CD2", "HD23" },
		{ "HG", "CG" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CE", "CD" },
		{ "CD", "CG" },
		{ "CD", "HD2" },
		{ "CD", "HD3" },
		{ "CE", "HE2" },
		{ "CE", "HE3" },
		{ "NZ", "CE" },
		{ "HG2", "CG" },
		{ "HG3", "CG" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
		{ "HZ1", "NZ" },
		{ "HZ2", "NZ" },
		{ "HZ3", "NZ" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "HE1", "CE" },
		{ "HE2", "CE" },
		{ "HE3", "CE" },
		{ "SD", "CG" },
		{ "CE", "SD" },
		{ "CG", "HG2" },
		{ "CG", "HG3" },
		{ "CG", "SD" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CE1", "CD1" },
		{ "CD1", "CG" },
		{ "HD1", "CD1" },
		{ "CE2", "CD2" },
		{ "CD2", "CG" },
		{ "HD2", "CD2" },
		{ "CE1", "CZ" },
		{ "CZ", "CE2" },
		{ "HE2", "CE2" },
		{ "HZ", "CZ" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD", "CG" },
		{ "CD", "HD2" },
		{ "CD", "HD3" },
		{ "CD", "N" },
		{ "CG", "HG2" },
		{ "CG", "HG3" },
		{ "H", "N" },
		{ "H3", "N" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "OG", "CB" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HG", "OG" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG2", "CB" },
		{ "HB", "CB" },
		{ "OG1", "CB" },
		{ "HG21", "CG2" },
		{ "HG22", "CG2" },
		{ "HG23", "CG2" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HG1", "OG1" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CD1", "CG" },
		{ "CD1", "HD1" },
		{ "NE1", "CD1" },
		{ "CE3", "CD2" },
		{ "CD2", "CG" },
		{ "CZ2", "CE2" },
		{ "CE2", "NE1" },
		{ "CZ3", "CE3" },
		{ "HE3", "CE3" },
		{ "CH2", "CZ2" },
		{ "CH2", "CZ3" },
		{ "HH2", "CH2" },
		{ "HZ2", "CZ2" },
		{ "HZ3", "CZ3" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HE1", "NE1" },
		{ "HXT", "OXT" },
	},
	{
		{ "N", "-C" },
		{ "C", "CA" },
		{ "O", "C" },
		{ "OXT", "C" },
		{ "CB", "CA" },
		{ "HA", "CA" },
		{ "CA", "N" },
		{ "CG", "CB" },
		{ "HB2", "CB" },
		{ "HB3", "CB" },
		{ "CE1", "CD1" },
		{ "CD1", "CG" },
		{ "HD1", "CD1" },
		{ "CE2", "CD2" },
		{ "CD2", "CG" },
		{ "CD2", "HD2" },
		{ "CE1", "CD1" },
		{ "CZ", "CE1" },
		{ "HE1", "CE1" },
		{ "CE2", "CD2" },
		{ "HE2", "CE2" },
		{ "OH", "CZ" },
		{ "H", "N" },
		{ "H2", "N" },
		{ "H3", "N" },
		{ "HH", "OH" },
		{ "HXT", "OXT" },
	},
	{ { "N", "-C" }, { "C", "CA" }, { "O", "C" }, { "OXT", "C" }, { "CB", "CA" }, { "HA", "CA" }, { "CA", "N" }, { "CG1", "CB" }, { "CG2", "CB" }, { "HB", "CB" }, { "CG1", "HG11" }, { "CG1", "HG12" }, { "CG1", "HG13" }, { "CG2", "HG21" }, { "CG2", "HG22" }, { "CG2", "HG23" }, { "H", "N" }, { "H2", "N" }, { "H3", "N" }, { "HXT", "OXT" } }
};

const static stringIndexMap aminoAcidIndex {
	{ "ARG", 0 },
	{ "ALA", 1 },
	{ "ASN", 2 },
	{ "ASP", 3 },
	{ "CYS", 4 },
	{ "GLN", 5 },
	{ "GLU", 6 },
	{ "GLY", 7 },
	{ "HIS", 8 },
	{ "ILE", 9 },
	{ "LEU", 10 },
	{ "LYS", 11 },
	{ "MET", 12 },
	{ "PHE", 13 },
	{ "PRO", 14 },
	{ "SER", 15 },
	{ "THR", 16 },
	{ "TRP", 17 },
	{ "TYR", 18 },
	{ "VAL", 19 }
};

/**
 *  The Residue class represent one of the twenty standard amino acids.
 *  A Residue is made of Atom objects which are connected via Bond objects.
 *
 *  @param[in] atoms_ Vector of Atom objects
 *  @param[in] aminoAcidName_ name of amino acid
 *  @param[in] residueName_ name of residue
 */
template <typename T>
Residue<T>::Residue(atomVector<T> atoms_, const std::string& aminoAcidName_, const std::string& residueName_, bool n_terminus, bool c_terminus)
{
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

	m_n_terminus = n_terminus;
	m_c_terminus = c_terminus;
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

template class prostruct::Residue<float>;
template class prostruct::Residue<double>;
