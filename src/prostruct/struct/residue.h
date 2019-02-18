//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_RESIDUE_H
#define PROSTRUCT_RESIDUE_H

#include "prostruct/struct/atom.h"
#include <armadillo>

using aminoAcidAtomMap = std::vector<std::map<std::string, std::string>>;
using stringIndexMap = std::map<std::string, int>;

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

template <typename T>
class Residue {
public:
	// takes an arbitrary number of atoms and tries to form a residue
	//    Residue(std::unique_ptr<Atom> atoms...);
	Residue(atomVector<T>, const std::string&, const std::string&);

	atomVector<T> getBackbone()
	{
		std::vector<std::shared_ptr<Atom<T>>> backboneAtoms;
		for (auto const& i : backbone) {
			backboneAtoms.emplace_back(atoms[i]);
		}
		return backboneAtoms;
	}
	atomVector<T> getSidechain()
	{
		std::vector<std::shared_ptr<Atom<T>>> sidechainAtoms;
		for (auto const& i : sidechain) {
			sidechainAtoms.emplace_back(atoms[i]);
		}
		return sidechainAtoms;
	}

	std::shared_ptr<Atom<T>> operator[](const int index) { return atoms[index]; }

	arma::Mat<T> getXYZ() { return xyz; }

	std::shared_ptr<Atom<T>> getAtom(int index) { return atoms[0]; }

	std::string getResidueName() { return residueName; }

	void link(std::shared_ptr<Residue<T>>);

	void createBonds();

	int n_atoms() { return atoms.size(); };

	void print_atoms()
	{

		for (auto const& atom : atoms) {
			std::cout << atom->getName() << "-";
		}
	}

	atomVector<T> getAtoms() { return atoms; }
	arma::Col<T> getRadii() { return radii; }

private:
	arma::Mat<T> xyz;
	arma::Col<T> radii;
	std::vector<int> backbone; /**< A vector with the index number of the backbone atoms */
	std::vector<int> sidechain; /**< A vector with the index number of the sidechain atoms */
	std::string aminoAcidName; /**< Name of the amino acid, e.g. ALA */
	enum AminoAcid aminoAcid; /**< Amino acid enum, e.g. ALA */
	std::string residueName; /**< Name of the residue, e.g. ALA1 */
	atomVector<T> atoms; /**< A vector with the pointers to the Atom objects */
	std::map<std::string, int> atomMap; /**< Map atom name to internal index */
};

template <typename T>
using residueVector = std::vector<std::shared_ptr<Residue<T>>>;

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

#endif //PROSTRUCT_RESIDUE_H
