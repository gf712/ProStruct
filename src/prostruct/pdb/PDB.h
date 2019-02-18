//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#ifndef NDEBUG
#define FMT_STRING_ALIAS 1
#endif

#include "prostruct/struct/chain.h"
#include "prostruct/parsers/PDBparser.h"
#include "prostruct/pdb/geometry.h"
#include "prostruct/utils/io.h"

#include <fmt/format.h>

namespace prostruct {
	template<typename T>
	class PDB {

	public:

		PDB(const std::string &filename);

		~PDB() = default;

		static PDB fetch(std::string);

        std::string to_string() {
            return format(fmt("<prostruct.PDB {} precision, with {} atoms, {} residues at {}>"),
                          demangled_type<T>(), m_natoms, m_nresidues, fmt::ptr(this));
        }

		int n_chains() { return numberOfChains; }

		std::vector<std::string> getChainIDs() {
			return chainOrder;
		}

		std::shared_ptr<Chain<T>> getChain(std::string name_) { return chainMap[name_]; }

		arma::Mat<T> calculate_KabschSander();

		arma::Mat<T> predict_backboneHbonds();

		void calculate_dssp();

		arma::Mat<T> getXYZ() { return xyz; }

		int n_residues() { return nResidues; }

		int n_atoms() { return nAtoms; }

		arma::Col<T> getRadii() { return radii; }

		arma::Col<T> calculate_ASA(T probe);

		double calculate_RMSD(PDB &other);

		arma::Col<T> calculate_centroid();
//    arma::Mat<T> select(std::string);

		void recentre();

		arma::Mat<T> calculate_phi_psi();

		void kabsch_rotation(PDB<T> &other);

		T kabsch_rmsd(PDB<T> &other);

//    void rotate(arma::Col<T> &rotation); // rotation = [rotation_x, rotation_y, rotation_z]
//    void rotate(T rotation_angle, std::string axis); // axis = {"x", "y", "z"}

	private:

		arma::Mat<T> xyz;
		std::string filename;
		int numberOfChains;
		std::map<std::string, std::shared_ptr<Chain<T>>> chainMap;
		int nAtoms;
		std::vector<std::string> chainOrder;
		arma::uword nResidues;
		arma::Col<T> radii;

		void internalKS(arma::Mat<T> &);
		void getBackboneAtoms(arma::Mat<T> &, arma::Mat<T> &, arma::Mat<T> &, arma::Mat<T> &);
	};

}

#endif //PROSTRUCT_PDB_H
