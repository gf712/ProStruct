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
	std::set<std::string> backbone_atom_names = { "C", "CA", "N", "O" };
	template <typename T> class PDB {
	public:
		PDB(const std::string& filename);

		static PDB fetch(std::string);

		std::string to_string()
		{
			return format(fmt("<prostruct.PDB {} precision, with {} atoms, {} "
							  "residues at {}>"),
				demangled_type<T>(), m_natoms, m_nresidues, fmt::ptr(this));
		}

		std::string get_filename() { return m_filename; }

		int n_chains() { return m_number_of_chains; }

		std::vector<std::string> get_chain_names() const noexcept
		{
			return m_chain_order;
		}

		std::shared_ptr<Chain<T>> get_chain(const std::string& name) const
		{
			return m_chain_map.at(name);
		}

		arma::Mat<T> compute_kabsch_sander() const;

		arma::Mat<T> predict_backbone_hbonds() const;

		void compute_dssp() const;

		arma::Mat<T> get_xyz() const noexcept { return m_xyz; }

		int n_residues() const noexcept { return m_nresidues; }

		int n_atoms() const noexcept { return m_natoms; }

		arma::Col<T> get_radii() const noexcept { return m_radii; }

		arma::Col<T> compute_shrake_rupley(
			T probe = 1.4, int n_sphere_points = 960) const;

		T calculate_RMSD(PDB& other) const;

		arma::Col<T> calculate_centroid() const;

		void recentre();

		arma::Mat<T> calculate_phi_psi(bool use_radians = false) const;

		void kabsch_rotation(PDB<T>& other);

		T kabsch_rmsd(PDB<T>& other) const;

		arma::Col<T> calculate_phi(bool use_radians = false) const;
		arma::Col<T> calculate_psi(bool use_radians = false) const;

		//    void rotate(arma::Col<T> &rotation); // rotation = [rotation_x,
		//    rotation_y, rotation_z] void rotate(T rotation_angle, std::string
		//    axis);
		//    // axis = {"x", "y", "z"}
#ifndef SWIG
		template <typename... Args>
		arma::Col<arma::uword> get_atom_indices(Args... patterns)
		{
			arma::Col<arma::uword> max_result(m_nresidues * m_natoms);
			arma::uword i = 0;
			arma::uword position_in_pdb = 0;

			for (const auto& chain_name : m_chain_order) {
				arma::Col<arma::uword> residue_atoms
					= m_chain_map[chain_name]->get_atom_indices(
						std::forward<Args>(patterns)...);
				max_result(arma::span(i, i + residue_atoms.n_rows - 1))
					= residue_atoms + position_in_pdb;
				i += residue_atoms.n_rows;
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
		residueVector<T> m_residues;
		static constexpr T to_rad_constant = 180.0 / M_PI;

		void internalKS(arma::Mat<T>&) const;
		arma::Mat<T> get_backbone_atoms() const;

		void append_new_residue(
			const std::pair<const std::string, atomVector<T>>&,
			residueVector<T>&, bool, bool);
	};

} // namespace prostruct

#endif // PROSTRUCT_PDB_H
