/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_KERNELS_H
#define PROSTRUCT_KERNELS_H

#include <prostruct/struct/utils.h>

#include <armadillo>

namespace prostruct
{
	// forward declare residue so don't have to add the header file here
	template <typename T>
	class Residue;
}

namespace prostruct::kernels
{
	/**
	 * The factor to convert radians to degrees
	 *
	 * @tparam T the required return type
	 */
	template <typename T>
	static constexpr T to_rad_constant = 180.0 / M_PI;
	/**
	 * An untyped version of dihedrals for armadillo optimisations
	 *
	 * @tparam T1 the required return type
	 * @tparam T2 the rvalue reference type of the atoms
	 * @param atom1 atom in the first plane
	 * @param atom2 atom in the first plane
	 * @param atom3 atom in the first plane and second plane
	 * @param atom4 atom in the second plane
	 * @param coef the coefficient to convert rad to degrees
	 * @return the dihedral angle
	 */
	template <typename T1, typename T2>
	inline T2 dihedrals_lazy(T1&& atom1, T1&& atom2, T1&& atom3, T1&& atom4, T2 coef)
	{
		arma::Col<T2> b1 = arma::normalise(atom1 - atom2);
		arma::Col<T2> b2 = arma::normalise(atom2 - atom3);
		arma::Col<T2> b3 = arma::normalise(atom3 - atom4);
		arma::Col<T2> n1 = arma::cross(b1, b2);
		arma::Col<T2> n2 = arma::cross(b2, b3);
		return std::atan2(arma::dot(arma::cross(n1, b2), n2), arma::dot(n1, n2)) * coef;
	}
	/**
	 * An untyped version of the norm for armadillo
	 * optimisations.
	 *
	 * @tparam T1 the required return type
	 * @tparam T2
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	template <typename T1, typename T2>
	inline T1 norm_lazy(T2&& atom1, T2&& atom2)
	{
		return arma::dot(atom1 - atom2, atom1 - atom2);
	}

	/**
	 * An untyped version of euclidean distances for armadillo
	 * optimisations.
	 * Note that the dot product tends to be quicker
	 * than arma::norm when using the MKL backend.
	 */
	template <typename T1, typename T2>
	inline T1 distance_lazy(T2&& atom1, T2&& atom2)
	{
		return std::sqrt(norm_lazy<T1>(atom1, atom2));
	}
	/**
	 * The phi dihedral angle kernel factory to calculate the
	 * phi angle formed by two residues.
	 *
	 * @tparam T the required return type
	 * @param use_radians whether to use radians or degrees
	 * @return lambda to perform phi_kernel
	 */
	template <typename T>
	auto phi_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto phi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
							  const std::shared_ptr<Residue<T>>& residue_next) -> T {
			if (residue->is_c_terminus())
				return 0.0;
			auto atom_coords_this = residue->get_backbone_atoms();
			auto atom_coords_next = residue_next->get_backbone_atoms();
			return kernels::dihedrals_lazy(atom_coords_this.col(2), atom_coords_next.col(0),
				atom_coords_next.col(1), atom_coords_next.col(2), coef);
		};
		return phi_kernel;
	}

	template <typename T>
	auto psi_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto psi_kernel = [coef](const std::shared_ptr<Residue<T>>& residue,
							  const std::shared_ptr<Residue<T>>& residue_next) -> T {
			if (residue_next->is_n_terminus())
				return 0.0;
			auto atom_coords_this = residue->get_backbone_atoms();
			auto atom_coords_next = residue_next->get_backbone_atoms();
			return kernels::dihedrals_lazy(atom_coords_this.col(0), atom_coords_this.col(1),
				atom_coords_this.col(2), atom_coords_next.col(0), coef);
		};
		return psi_kernel;
	}
	template <typename T>
	auto chi1_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto chi1_kernel = [coef](const std::shared_ptr<Residue<T>>& residue) -> T {
			arma::Mat<T> coords;
			switch (residue->get_amino_acid_type())
			{
			case AminoAcid::VAL:
			case AminoAcid::ILE:
			{
				coords = residue->get_atom_coords("N", "CA", "CB", "CG1");
			}
			break;
			case AminoAcid::THR:
			{
				coords = residue->get_atom_coords("N", "CA", "CB", "OG1");
			}
			break;
			case AminoAcid::SER:
			{
				coords = residue->get_atom_coords("N", "CA", "CB", "OG");
			}
			break;
			case AminoAcid::CYS:
			{
				coords = residue->get_atom_coords("N", "CA", "CB", "SG");
			}
			break;
			case AminoAcid::ALA:
			case AminoAcid::GLY:
			{
				return 0.0;
			}
			default:
			{
				coords = residue->get_atom_coords("N", "CA", "CB", "CG");
			}
			}

			return kernels::dihedrals_lazy(
				coords.col(0), coords.col(1), coords.col(2), coords.col(3), coef);
		};
		return chi1_kernel;
	}

	template <typename T>
	auto chi2_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		// the chi2 kernel as a C++ lambda
		auto chi2_kernel = [coef](const std::shared_ptr<Residue<T>>& residue) -> T {
			arma::Mat<T> coords;
			switch (residue->get_amino_acid_type())
			{
			case AminoAcid::ARG:
			case AminoAcid::GLN:
			case AminoAcid::GLU:
			case AminoAcid::LYS:
			case AminoAcid::PRO:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG", "CD");
			}
			break;
			case AminoAcid::ASN:
			case AminoAcid::ASP:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG", "OD1");
			}
			break;
			case AminoAcid::HIS:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG", "ND1");
			}
			break;
			case AminoAcid::ILE:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG1", "CD1");
			}
			break;
			case AminoAcid::LEU:
			case AminoAcid::PHE:
			case AminoAcid::TRP:
			case AminoAcid::TYR:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG", "CD1");
			}
			break;
			case AminoAcid::MET:
			{
				coords = residue->get_atom_coords("CA", "CB", "CG", "SD");
			}
			break;
			default:
			{
				return 0.0;
			}
			}

			return kernels::dihedrals_lazy(
				coords.col(0), coords.col(1), coords.col(2), coords.col(3), coef);
		};
		return chi2_kernel;
	}
	template <typename T>
	auto chi3_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto chi3_kernel = [coef](const std::shared_ptr<Residue<T>>& residue) -> T {
			arma::Mat<T> coords;
			switch (residue->get_amino_acid_type())
			{
			case AminoAcid::ARG:
			{
				coords = residue->get_atom_coords("CB", "CG", "CD", "NE");
			}
			break;
			case AminoAcid::GLN:
			case AminoAcid::GLU:
			{
				coords = residue->get_atom_coords("CB", "CG", "CD", "OE1");
			}
			break;
			case AminoAcid::LYS:
			{
				coords = residue->get_atom_coords("CB", "CG", "CD", "CE");
			}
			break;
			case AminoAcid::MET:
			{
				coords = residue->get_atom_coords("CB", "CG", "SD", "CE");
			}
			break;
			default:
			{
				return 0.0;
			}
			}
			return kernels::dihedrals_lazy(
				coords.col(0), coords.col(1), coords.col(2), coords.col(3), coef);
		};
		return chi3_kernel;
	}

	template <typename T>
	auto chi4_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto chi4_kernel = [coef](const std::shared_ptr<Residue<T>>& residue) -> T {
			arma::Mat<T> coords;
			switch (residue->get_amino_acid_type())
			{
			case AminoAcid::ARG:
			{
				coords = residue->get_atom_coords("CG", "CD", "NE", "CZ");
			}
			break;
			case AminoAcid::LYS:
			{
				coords = residue->get_atom_coords("CG", "CD", "CE", "NZ");
			}
			break;
			default:
			{
				return 0.0;
			}
			}
			return kernels::dihedrals_lazy(
				coords.col(0), coords.col(1), coords.col(2), coords.col(3), coef);
		};
		return chi4_kernel;
	}

	template <typename T>
	auto chi5_kernel(bool use_radians)
	{
		T coef = use_radians ? 1.0 : to_rad_constant<T>;
		auto chi5_kernel = [coef](const std::shared_ptr<Residue<T>>& residue) -> T {
			arma::Mat<T> coords;
			switch (residue->get_amino_acid_type())
			{
			case AminoAcid::ARG:
			{
				coords = residue->get_atom_coords("CD", "NE", "CZ", "NH1");
			}
			break;
			default:
			{
				return 0.0;
			}
			}
			return kernels::dihedrals_lazy(
				coords.col(0), coords.col(1), coords.col(2), coords.col(3), coef);
		};
		return chi5_kernel;
	}
}

#endif // PROSTRUCT_KERNELS_H
