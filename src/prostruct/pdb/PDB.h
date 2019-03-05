//
// Created by gil on 27/03/18.
//

#ifndef PROSTRUCT_PDB_H
#define PROSTRUCT_PDB_H

#define FMT_STRING_ALIAS 1

#include "prostruct/parsers/PDBparser.h"
#include "prostruct/pdb/geometry.h"
#include "prostruct/struct/chain.h"
#include "prostruct/utils/io.h"

#include <set>

#include <fmt/format.h>

namespace prostruct {
    std::set<std::string> backbone_atom_names = {"C", "CA", "N", "O"};
template <typename T>
class PDB {
public:
	PDB(const std::string& filename);

	~PDB() = default;

	static PDB fetch(std::string);

	std::string to_string()
	{
		return format(
			fmt("<prostruct.PDB {} precision, with {} atoms, {} residues at {}>"),
			demangled_type<T>(), m_natoms, m_nresidues, fmt::ptr(this));
	}

	std::string get_filename() { return m_filename; }

	int n_chains() { return m_number_of_chains; }

	std::vector<std::string> get_chain_names() { return m_chain_order; }

	std::shared_ptr<Chain<T>> get_chain(std::string name_)
	{
		return m_chain_map[name_];
	}

	arma::Mat<T> compute_kabsch_sander();

	arma::Mat<T> predict_backbone_hbonds();

	void compute_dssp();

	arma::Mat<T> get_xyz() { return m_xyz; }

	int n_residues() { return m_nresidues; }

	int n_atoms() { return m_natoms; }

	arma::Col<T> get_radii() { return m_radii; }

	arma::Col<T> compute_shrake_rupley(T probe=1.4, int n_sphere_points=960);

	T calculate_RMSD(PDB& other);

	arma::Col<T> calculate_centroid();

	void recentre();

	arma::Mat<T> calculate_phi_psi();

	void kabsch_rotation(PDB<T>& other);

	T kabsch_rmsd(PDB<T>& other);

	//    void rotate(arma::Col<T> &rotation); // rotation = [rotation_x,
	//    rotation_y, rotation_z] void rotate(T rotation_angle, std::string axis);
	//    // axis = {"x", "y", "z"}
#ifndef SWIG
	template <typename...Args>
	arma::Col<arma::uword> get_atom_indices(Args... patterns) 
	{
		arma::Col<arma::uword> max_result(m_nresidues*m_natoms);
		arma::uword i = 0;
	    arma::uword position_in_pdb = 0;
	    
		for(const auto& chain_name: m_chain_order) {
			arma::Col<arma::uword> residue_atoms = m_chain_map[chain_name]->get_atom_indices(std::forward<Args>(patterns)...);
			max_result(arma::span(i, i+residue_atoms.n_rows-1)) = residue_atoms + position_in_pdb;
			i+=residue_atoms.n_rows;
			position_in_pdb += m_chain_map[chain_name]->m_natoms;
		}

		arma::Col<arma::uword> result = max_result.head(i);

		return result;
	}
#endif

private:
	arma::Mat<T> m_xyz;
	std::string m_filename;
	int m_number_of_chains;
	std::map<std::string, std::shared_ptr<Chain<T>>> m_chain_map;
	int m_natoms;
	std::vector<std::string> m_chain_order;
	arma::uword m_nresidues;
	arma::Col<T> m_radii;

	void internalKS(arma::Mat<T>&);
	arma::Mat<T> get_backbone_atoms();
};

} // namespace prostruct

#endif // PROSTRUCT_PDB_H
