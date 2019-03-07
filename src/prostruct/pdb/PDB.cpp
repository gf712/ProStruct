/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include <prostruct/core/engine.h>
#include <prostruct/core/kernels.h>
#include <prostruct/pdb/PDB.h>

using namespace prostruct;

template <typename T> PDB<T>::PDB(const std::string& filename_)
{
	m_filename = filename_;

	std::map<std::string, std::map<std::string, atomVector<T>, AASequenceOrder>>
		chainAtomMap;

	createMap(m_filename, chainAtomMap, m_chain_order);

	m_natoms = 0;
	m_nresidues = 0;
	m_xyz.set_size(3, 0);

	for (const auto& chain : m_chain_order) {
		residueVector<T> residues;
		const std::map<std::string, atomVector<T>, AASequenceOrder> chain_i
			= chainAtomMap.at(chain);
		append_new_residue(*chain_i.cbegin(), residues, true, false);
		for (auto atomPair = std::next(chain_i.cbegin());
			 atomPair != std::prev(chain_i.cend()); ++atomPair) {
			append_new_residue(*atomPair, residues, false, false);
		}
		append_new_residue(*std::prev(chain_i.cend()), residues, false, true);
		m_chain_map[chain] = std::make_shared<Chain<T>>(residues, chain);
		m_nresidues += static_cast<arma::uword>(residues.size());
		m_residues.insert(m_residues.end(), residues.begin(), residues.end());
	}
	m_number_of_chains = static_cast<int>(m_chain_map.size());
}

template <typename T>
void PDB<T>::append_new_residue(
	const std::pair<const std::string, atomVector<T>>& atom_pair,
	residueVector<T>& residues, bool n_terminus, bool c_terminus)
{
	residues.emplace_back(std::make_shared<Residue<T>>(atom_pair.second,
		atom_pair.first.substr(0, 3), atom_pair.first, n_terminus, c_terminus));
	m_xyz.insert_cols(
		static_cast<arma::uword>(m_natoms), residues.back()->getXYZ());
	m_radii.insert_rows(
		static_cast<arma::uword>(m_natoms), residues.back()->getRadii());
	m_natoms += residues.back()->n_atoms();
}

template <typename T> PDB<T> PDB<T>::fetch(std::string PDB_id)
{
	throw "Not implemented";
}

template <typename T> arma::Mat<T> PDB<T>::get_backbone_atoms() const
{
	arma::Mat<T> result(3, m_nresidues * 4);
	arma::uword pos = 0;
	arma::uword i = 0;

	for (auto const& chain : m_chain_order) {
		for (auto const& residue : m_chain_map.at(chain)->getResidues()) {
			result(arma::span::all, arma::span(i, i + 3))
				= m_xyz(arma::span::all, arma::span(pos, pos + 3));
			pos += residue->n_atoms();
			i += 4;
		}
	}
	return result;
}

template <typename T> void PDB<T>::internalKS(arma::Mat<T>& E) const
{
	auto backbone_atom_coords = get_backbone_atoms();
	geometry::kabsch_sander(backbone_atom_coords, E);
}

template <typename T> arma::Mat<T> PDB<T>::compute_kabsch_sander() const
{
	arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);

	internalKS(E);

	return E;
};

template <typename T>
arma::Col<T> PDB<T>::compute_shrake_rupley(T probe, int n_sphere_points) const
{
	// calculates atom surface accessibility using the Shrake-Rupley algorithm

	arma::Col<T> asa(static_cast<const arma::uword>(m_natoms));

	geometry::shrake_rupley(m_xyz, m_radii, asa,
		static_cast<arma::uword>(m_natoms), probe,
		static_cast<arma::uword>(n_sphere_points));

	return asa;
}

template <typename T> arma::Mat<T> PDB<T>::predict_backbone_hbonds() const
{
	arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);

	internalKS(E);

	E.for_each([](T& elem) { elem = elem < -0.5; });

	return E;
}

template <typename T> void PDB<T>::compute_dssp() const
{
	// arma::Mat<T> C_coords(3, m_nresidues);
	// arma::Mat<T> O_coords(3, m_nresidues);
	// arma::Mat<T> N_coords(3, m_nresidues);
	// arma::Mat<T> CA_coords(3, m_nresidues);

	// get_backbone_atoms(C_coords, O_coords, N_coords, CA_coords);

	auto backbone_atoms = get_backbone_atoms();

	// geometry::dssp(C_coords, O_coords, N_coords, CA_coords);
	// geometry::dssp(xyz, backbone_atoms);
}

template <typename T> T PDB<T>::calculate_RMSD(PDB& other) const
{
	// first check if the size is the same
	if (m_natoms != other.n_atoms()) {
		throw "Atom number mismatch";
	}

	return geometry::rmsd(m_xyz, other.get_xyz());
}

template <typename T> arma::Col<T> PDB<T>::calculate_centroid() const
{
	arma::Col<T> result(3);

	geometry::get_centroid(m_xyz, result);

	return result;
}

template <typename T> void PDB<T>::recentre()
{
	geometry::recentre_molecule(m_xyz);
}

template <typename T>
arma::Mat<T> PDB<T>::calculate_phi_psi(bool use_radians) const
{
	T coef = use_radians ? 1.0 : to_rad_constant;

	// the Phi kernel as a C++ lambda
	auto phi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
						  const std::shared_ptr<Residue<T>>& residue_next) {
		if (residue->is_c_terminus())
			return static_cast<T>(0.0);
		auto atom_coords_this = residue->get_backbone_atoms();
		auto atom_coords_next = residue_next->get_backbone_atoms();
		return kernels::dihedrals_lazy(atom_coords_this.col(2),
			atom_coords_next.col(0), atom_coords_next.col(1),
			atom_coords_next.col(2), coef);
	};

	// the Psi kernel as a C++ lambda
	auto psi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
						  const std::shared_ptr<Residue<T>>& residue_next) {
		if (residue_next->is_n_terminus())
			return static_cast<T>(0.0);
		auto atom_coords_this = residue->get_backbone_atoms();
		auto atom_coords_next = residue_next->get_backbone_atoms();
		return kernels::dihedrals_lazy(atom_coords_this.col(0),
			atom_coords_this.col(1), atom_coords_this.col(2),
			atom_coords_next.col(0), coef);
	};

	return core::residue_kernel_engine(m_residues, 0, phi_kernel, psi_kernel);
}

template <typename T> arma::Col<T> PDB<T>::calculate_phi(bool use_radians) const
{
	// convert from radians to degrees
	T coef = use_radians ? 1.0 : to_rad_constant;
	// the Phi kernel as a C++ lambda
	auto phi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
						  const std::shared_ptr<Residue<T>>& residue_next) {
		if (residue->is_c_terminus())
			return static_cast<T>(0.0);
		auto atom_coords_this = residue->get_backbone_atoms();
		auto atom_coords_next = residue_next->get_backbone_atoms();
		return kernels::dihedrals_lazy(atom_coords_this.col(2),
			atom_coords_next.col(0), atom_coords_next.col(1),
			atom_coords_next.col(2), coef);
	};

	// result is a matrix where each row has the lambda/kernel result for each
	// residue, so transform it into column vector
	return arma::Col<T>(
		core::residue_kernel_engine(m_residues, 0, phi_kernel).memptr(),
		m_nresidues);
}

template <typename T> arma::Col<T> PDB<T>::calculate_psi(bool use_radians) const
{
	// convert from radians to degrees
	T coef = use_radians ? 1.0 : to_rad_constant;
	// the Psi kernel as a C++ lambda
	auto psi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
						  const std::shared_ptr<Residue<T>>& residue_next) {
		if (residue_next->is_n_terminus())
			return static_cast<T>(0.0);
		auto atom_coords_this = residue->get_backbone_atoms();
		auto atom_coords_next = residue_next->get_backbone_atoms();
		return kernels::dihedrals_lazy(atom_coords_this.col(0),
			atom_coords_this.col(1), atom_coords_this.col(2),
			atom_coords_next.col(0), coef);
	};

	// result is a matrix where each row has the lambda/kernel result for each
	// residue, so transform it into column vector
	return arma::Col<T>(
		core::residue_kernel_engine(m_residues, 0, psi_kernel).memptr(),
		m_nresidues);
}
// void rotate(arma::Col<T> &rotation) {
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

template <typename T> void PDB<T>::kabsch_rotation(PDB<T>& other)
{
	// make copy of xyz
	auto xyz_copy = other.get_xyz();
	geometry::kabsch_rotation_(m_xyz, xyz_copy);
}

template <typename T> T PDB<T>::kabsch_rmsd(PDB& other) const
{
	// make copy of xyz
	auto xyz_copy = m_xyz;
	auto xyz_other_copy = other.get_xyz();

	return geometry::kabsch_rmsd_(xyz_copy, xyz_other_copy);
}

template class prostruct::PDB<float>;
template class prostruct::PDB<double>;