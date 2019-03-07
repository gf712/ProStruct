//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_RESIDUE_H
#define PROSTRUCT_RESIDUE_H

#include <regex>

#include "prostruct/struct/atom.h"

namespace prostruct {

	using aminoAcidAtomMap = std::vector<std::map<std::string, std::string>>;
	using stringIndexMap = std::map<std::string, int>;

	enum class aaLocation { Backbone, Sidechain };

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

	template <typename T> class Residue {
		template <typename T1> friend class Chain;

	public:
		// takes an arbitrary number of atoms and tries to form a residue
		//    Residue(std::unique_ptr<Atom> atoms...);
		Residue(atomVector<T>, const std::string&, const std::string&,
			bool = false, bool = false);

		atomVector<T> getBackbone()
		{
			std::vector<std::shared_ptr<Atom<T>>> backboneAtoms;
			for (auto const& i : backbone) {
				backboneAtoms.emplace_back(atoms[i]);
			}
			return backboneAtoms;
		}

		arma::Mat<T> get_backbone_atoms()
		{
			return xyz(arma::span::all, arma::span(0, 3));
		}

		atomVector<T> getSidechain()
		{
			std::vector<std::shared_ptr<Atom<T>>> sidechainAtoms;
			for (auto const& i : sidechain) {
				sidechainAtoms.emplace_back(atoms[i]);
			}
			return sidechainAtoms;
		}

		std::shared_ptr<Atom<T>> operator[](const int index)
		{
			return atoms[index];
		}

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

		bool is_n_terminus() { return m_n_terminus; }

		bool is_c_terminus() { return m_c_terminus; }

	protected:
		template <typename... Args>
		arma::Col<arma::uword> get_atom_indices(Args... patterns)
		{
			arma::Col<arma::uword> max_result(atoms.size());
			std::vector<std::string> vec = { patterns... };

			arma::uword result_idx = 0;
			for (const auto& pattern : vec) {
				arma::uword i = 0;
				for (const auto& atom : atoms) {
					if (atom->getName() == pattern) {
						max_result(result_idx) = i;
						result_idx++;
					}
					i++;
				}
			}
			arma::Col<arma::uword> result = max_result.head(result_idx);
			return result;
		}

	private:
		arma::Mat<T> xyz;
		arma::Col<T> radii;
		bool m_n_terminus;
		bool m_c_terminus;
		std::vector<int> backbone; /**< A vector with the index number of the
									  backbone atoms */
		std::vector<int> sidechain; /**< A vector with the index number of the
									   sidechain atoms */
		std::string aminoAcidName; /**< Name of the amino acid, e.g. ALA */
		enum AminoAcid aminoAcid; /**< Amino acid enum, e.g. ALA */
		std::string residueName; /**< Name of the residue, e.g. ALA1 */
		atomVector<T>
			atoms; /**< A vector with the pointers to the Atom objects */
		std::map<std::string, int>
			atomMap; /**< Map atom name to internal index */
	};

	template <typename T>
	using residueVector = std::vector<std::shared_ptr<Residue<T>>>;
}

#endif // PROSTRUCT_RESIDUE_H
