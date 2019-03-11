/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#ifndef PROSTRUCT_STRUCT_BASE_H
#define PROSTRUCT_STRUCT_BASE_H

#define FMT_STRING_ALIAS 1

#include <prostruct/core/engine.h>
#include <prostruct/core/kernels.h>
#include <prostruct/pdb/geometry.h>
#include <prostruct/struct/residue.h>
#include <prostruct/utils/io.h>

#include <fmt/format.h>

using namespace prostruct;

namespace prostruct
{
	template <class T>
	class StructBase
	{
	public:
		StructBase() {};

		virtual std::string to_string() const
		{
			return format(fmt("<prostruct.StructBase {} precision, with {} atoms, {} "
							  "residues at {}>"),
				demangled_type<T>(), m_natoms, m_nresidues, fmt::ptr(this));
		}

		arma::Mat<T> compute_kabsch_sander() const noexcept
		{
			arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);
			internalKS(E);
			return E;
		}

		arma::Mat<T> predict_backbone_hbonds() const noexcept
		{
			arma::Mat<T> E(m_nresidues, m_nresidues, arma::fill::zeros);
			internalKS(E);
			E.for_each([](T& elem) { elem = elem < -0.5; });
			return E;
		}

		void compute_dssp() const noexcept
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

		arma::Mat<T> get_xyz() const noexcept { return m_xyz; }

		int n_residues() const noexcept { return m_nresidues; }

		int n_atoms() const noexcept { return m_natoms; }

		arma::Col<T> get_radii() const noexcept { return m_radii; }

		std::vector<std::shared_ptr<Residue<T>>> get_residues() const noexcept
		{
			return m_residues;
		}

		arma::Col<T> compute_shrake_rupley(T probe = 1.4, int n_sphere_points = 960) const noexcept
		{
			arma::Col<T> asa(static_cast<arma::uword>(m_natoms));

			geometry::shrake_rupley(m_xyz, m_radii, asa, static_cast<arma::uword>(m_natoms), probe,
				static_cast<arma::uword>(n_sphere_points));

			return asa;
		}

		T calculate_RMSD(StructBase<T>& other) const
		{
			// first check if the size is the same
			if (m_natoms != other.n_atoms())
			{
				throw "Atom number mismatch";
			}

			return geometry::rmsd(m_xyz, other.get_xyz());
		}

		arma::Col<T> calculate_centroid() const noexcept
		{
			arma::Col<T> result(3);
			geometry::get_centroid(m_xyz, result);
			return result;
		}

		void recentre() noexcept { geometry::recentre_molecule(m_xyz); }

		arma::Mat<T> calculate_phi_psi(bool use_radians = false) const noexcept
		{
			return core::residue_kernel_engine(m_residues, 0, kernels::phi_kernel<T>(use_radians),
				kernels::psi_kernel<T>(use_radians));
		}

		void kabsch_rotation(StructBase<T>& other) noexcept
		{
			// make copy of xyz
			auto xyz_copy = other.get_xyz();
			geometry::kabsch_rotation_(m_xyz, xyz_copy);
		}

		T kabsch_rmsd(StructBase<T>& other) const noexcept
		{
			// make copy of xyz
			auto xyz_copy = m_xyz;
			auto xyz_other_copy = other.get_xyz();

			return geometry::kabsch_rmsd_(xyz_copy, xyz_other_copy);
		}

		arma::Col<T> calculate_phi(bool use_radians = false) const noexcept
		{
			// result is a matrix where each row has the lambda/kernel result for each
			// residue, so transform it into column vector
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::phi_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_psi(bool use_radians = false) const noexcept
		{
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::psi_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_chi1(bool use_radians = false) const noexcept
		{
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::chi1_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_chi2(bool use_radians = false) const noexcept
		{
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::chi2_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_chi3(bool use_radians = false) const noexcept
		{
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::chi3_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_chi4(bool use_radians = false) const noexcept
		{
			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::chi4_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Col<T> calculate_chi5(bool use_radians = false) const noexcept
		{

			return arma::Col<T>(
				core::residue_kernel_engine(m_residues, 0, kernels::chi5_kernel<T>(use_radians))
					.memptr(),
				m_nresidues);
		}

		arma::Mat<T> compute_shortest_distance() const noexcept
		{
			auto neighbour_kernel = [](const std::shared_ptr<Residue<T>>& residue,
										const std::shared_ptr<Residue<T>>& residue_neighbour) {
				auto sidechain = residue->get_sidechain_atoms();
				auto sidechain_neighbour = residue_neighbour->get_sidechain_atoms();
				T shortest_distance = std::numeric_limits<T>::infinity();

				for (arma::uword i = 0; i < sidechain.n_cols; ++i)
				{
					for (arma::uword j = i + 1; j < sidechain_neighbour.n_cols; ++j)
					{
						auto dist = kernels::distance_lazy<T>(
							sidechain.col(i), sidechain_neighbour.col(j));
						if (dist < shortest_distance)
							shortest_distance = dist;
					}
				}
				return shortest_distance;
			};
			return core::pairwise_residue_kernel_engine(m_residues, 0, neighbour_kernel).slice(0);
		}

		arma::Mat<T> compute_neighbours(T threshold) const noexcept
		{
			auto neighbour_kernel = [threshold](const std::shared_ptr<Residue<T>>& residue,
										const std::shared_ptr<Residue<T>>& residue_neighbour) {
				auto sidechain = residue->get_sidechain_atoms();
				auto sidechain_neighbour = residue_neighbour->get_sidechain_atoms();
				T found_neighbour = 0;

				for (arma::uword i = 0; i < sidechain.n_cols; ++i)
				{
					for (arma::uword j = i + 1; j < sidechain_neighbour.n_cols; ++j)
					{
						auto dist = kernels::distance_lazy<T>(
							sidechain.col(i), sidechain_neighbour.col(j));
						if (dist < threshold)
							return 1;
					}
				}
				return 0;
			};
			return core::pairwise_residue_kernel_engine(m_residues, 0, neighbour_kernel).slice(0);
		}

		//    void rotate(arma::Col<T> &rotation); // rotation = [rotation_x,
		//    rotation_y, rotation_z] void rotate(T rotation_angle, std::string
		//    axis);
		//    // axis = {"x", "y", "z"}
		void append_new_residue(const std::pair<const std::string, atomVector<T>>& atom_pair,
			residueVector<T>& residues, bool n_terminus, bool c_terminus)
		{
			residues.emplace_back(std::make_shared<Residue<T>>(atom_pair.second,
				atom_pair.first.substr(0, 3), atom_pair.first, n_terminus, c_terminus));
			m_xyz.insert_cols(static_cast<arma::uword>(m_natoms), residues.back()->get_xyz());
			m_radii.insert_rows(static_cast<arma::uword>(m_natoms), residues.back()->getRadii());
			m_natoms += residues.back()->n_atoms();
		}

		void append_new_residue(const Residue<T>& residue, bool n_terminus, bool c_terminus)
		{
			m_xyz.insert_cols(static_cast<arma::uword>(m_natoms), residue.get_xyz());
			m_radii.insert_rows(static_cast<arma::uword>(m_natoms), residue.getRadii());
			m_natoms += residue.n_atoms();
		}

		arma::Mat<T> get_backbone_atoms() const noexcept
		{
			arma::Mat<T> result(3, m_nresidues * 4);
			arma::uword pos = 0;
			arma::uword i = 0;

			for (auto const& residue : m_residues)
			{
				result(arma::span::all, arma::span(i, i + 3))
					= m_xyz(arma::span::all, arma::span(pos, pos + 3));
				pos += residue->n_atoms();
				i += 4;
			}
			return result;
		}

	protected:
		arma::Mat<T> m_xyz;
		int m_natoms;
		arma::uword m_nresidues;
		arma::Col<T> m_radii;
		residueVector<T> m_residues;
		static constexpr T to_rad_constant = 180.0 / M_PI;
		void internalKS(arma::Mat<T>& E) const noexcept
		{
			auto backbone_atom_coords = get_backbone_atoms();
			geometry::kabsch_sander(backbone_atom_coords, E);
		}
	};
}

#endif // PROSTRUCT_STRUCT_BASE_H
