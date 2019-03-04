//
// Created by gil on 27/03/18.
//

#include "prostruct/pdb/PDB.h"
#include "PDB.h"


using namespace prostruct;

template <typename T>
PDB<T>::PDB(const std::string& filename_)
{
	m_filename = filename_;

	std::map<std::string, std::map<std::string, atomVector<T>, AASequenceOrder>> chainAtomMap;

	createMap(m_filename, chainAtomMap, m_chain_order);

	m_natoms = 0;
	m_nresidues = 0;

	m_xyz.resize(3, 0);

	for (const auto& chain : m_chain_order) {

		residueVector<T> residues;

		for (auto const& atomPair : chainAtomMap[chain]) {

			try {
				auto residue = std::make_shared<Residue<T>>(atomPair.second, atomPair.first.substr(0, 3), atomPair.first);
				residues.emplace_back(residue);
				m_xyz.insert_cols(static_cast<const arma::uword>(m_natoms), residue->getXYZ());
				m_radii.insert_rows(static_cast<const arma::uword>(m_natoms), residue->getRadii());
				m_natoms += residue->n_atoms();
			} catch (const char* msg) {
				std::cout << "Residue: " << atomPair.first << std::endl;
				std::cout << msg << std::endl;
			}
		}

		try {
			m_chain_map[chain] = std::make_shared<Chain<T>>(residues, chain);
			m_nresidues += static_cast<arma::uword>(residues.size());
		} catch (const char* msg) {
			std::cout << "Chain: " << chain << std::endl;
			std::cout << msg << std::endl;
		}
	}

	m_number_of_chains = static_cast<int>(m_chain_map.size());
}

template <typename T>
PDB<T> PDB<T>::fetch(std::string PDB_id)
{
	throw "Not implemented";
}

template <typename T>
arma::Mat<T> PDB<T>::get_backbone_atoms()
{
	arma::Mat<T> result(3, m_nresidues * 4);
	arma::uword pos = 0;
	arma::uword i = 0;

	for (auto const& chain : m_chain_order) {
		for (auto const& residue : m_chain_map[chain]->getResidues()) {
			result(arma::span::all, arma::span(i, i+3)) = m_xyz(arma::span::all, arma::span(pos, pos+3));
			pos+=residue->n_atoms();
			i+=4;
		}
	}
	return result;
}

template <typename T>
void PDB<T>::internalKS(arma::Mat<T>& E)
{
	auto backbone_atom_coords = get_backbone_atoms();
    geometry::kabsch_sander(backbone_atom_coords, E);
}

template <typename T>
arma::Mat<T> PDB<T>::compute_kabsch_sander()
{
	arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);

	internalKS(E);

	return E;
};

template <typename T>
arma::Col<T> PDB<T>::compute_shrake_rupley(T probe, int n_sphere_points)
{
	// calculates atom surface accessibility using the Shrake-Rupley algorithm

	arma::Col<T> asa(static_cast<const arma::uword>(m_natoms));

    geometry::shrake_rupley(m_xyz, m_radii, asa, static_cast<arma::uword>(m_natoms), probe,
            static_cast<arma::uword>(n_sphere_points));

	return asa;
}

template <typename T>
arma::Mat<T> PDB<T>::predict_backbone_hbonds()
{
	arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);

	internalKS(E);

	E.for_each([](T& elem) { elem = elem < -0.5; });

	return E;
}

template <typename T>
void PDB<T>::compute_dssp()
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

template <typename T>
T PDB<T>::calculate_RMSD(PDB& other)
{
	// first check if the size is the same
	if (m_natoms != other.n_atoms()) {
		throw "Atom number mismatch";
	}

	return geometry::rmsd(m_xyz, other.get_xyz());
}

template <typename T>
arma::Col<T> PDB<T>::calculate_centroid()
{
	arma::Col<T> result(3);

	geometry::get_centroid(m_xyz, result);

	return result;
}

template <typename T>
void PDB<T>::recentre()
{
    geometry::recentre_molecule(m_xyz);
}

template <typename T>
arma::Mat<T> PDB<T>::calculate_phi_psi()
{

	// calculate phi and psi dihedral angles along the protein backbone
	// returns a matrix of shape (2, n_atoms):
	//  - row 0: phi
	//  - row 1: psi

	// different chains are not connected, therefore need ignore the C and N terminus of consecutive chains

	arma::Mat<T> result(m_nresidues - m_chain_order.size() + 1, 2, arma::fill::zeros);

	//	arma::Mat<T> result(3, nResidues*4, arma::fill::zeros);

	arma::uword atomPosition = 0;

	auto backbone_atoms = get_backbone_atoms();

	for (const auto& chainName : m_chain_order) {
		auto chain = m_chain_map[chainName];
		arma::uword n_residues = chain->n_residues();

		// create a 3D tensor with all the relevant coordinates
		arma::Cube<T> phiAtomCoords(3, 4, n_residues - 1);
		arma::Cube<T> psiAtomCoords(3, 4, n_residues - 1);

		for (arma::uword i = atomPosition; i < n_residues - 1; ++i) {
			// C|O|N|CA
			// phi angle
			phiAtomCoords.slice(i)(arma::span::all, 0) = backbone_atoms.col(i*4);
			phiAtomCoords.slice(i)(arma::span::all, 1) = backbone_atoms.col((i + 1)*4+2);
			phiAtomCoords.slice(i)(arma::span::all, 2) = backbone_atoms.col((i + 1)*4+3);
			phiAtomCoords.slice(i)(arma::span::all, 3) = backbone_atoms.col((i + 1)*4);

			// psi angle
			psiAtomCoords.slice(i)(arma::span::all, 0) = backbone_atoms.col(i*4+2);
			psiAtomCoords.slice(i)(arma::span::all, 1) = backbone_atoms.col(i*4+3);
			psiAtomCoords.slice(i)(arma::span::all, 2) = backbone_atoms.col(i*4);
			psiAtomCoords.slice(i)(arma::span::all, 3) = backbone_atoms.col((i + 1)*4+2);
		}

		arma::Col<T> phi(n_residues - 1);
		arma::Col<T> psi(n_residues - 1);

        geometry::dihedrals(phiAtomCoords, phi);
        geometry::dihedrals(psiAtomCoords, psi);

		result(arma::span(atomPosition, n_residues - 2), 0) = phi * (180.0 / M_PI);
		result(arma::span(atomPosition + 1, n_residues - 1), 1) = psi * (180.0 / M_PI);

		atomPosition += n_residues + 1;
        result.print();
        return result;
	}
    result.print();
	return result;
}

//void rotate(arma::Col<T> &rotation) {
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

template <typename T>
void PDB<T>::kabsch_rotation(PDB<T>& other)
{
	// make copy of xyz
	auto xyz_copy = other.get_xyz();
    geometry::kabsch_rotation_(m_xyz, xyz_copy);
}

template <typename T>
T PDB<T>::kabsch_rmsd(PDB& other)
{
	// make copy of xyz
	auto xyz_copy = m_xyz;
	auto xyz_other_copy = other.get_xyz();

	return geometry::kabsch_rmsd_(xyz_copy, xyz_other_copy);
}

template class prostruct::PDB<float>;
template class prostruct::PDB<double>;