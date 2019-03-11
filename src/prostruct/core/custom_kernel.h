/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include <prostruct/core/engine.h>
#include <prostruct/pdb/PDB.h>
#include <prostruct/utils/io.h>

#ifndef PROSTRUCT_CustomPDB_H
#define PROSTRUCT_CustomPDB_H

#ifdef SWIGPYTHON

namespace prostruct
{
	class CustomPDB : public PDB<float>
	{
	public:
		CustomPDB(const std::string& filename): PDB<float>(filename) {};

		virtual ~CustomPDB() {}

		virtual std::string to_string() const noexcept
		{
			// repeated here due to swig director (bug?)
			return format(fmt("<prostruct.PDB {} precision, with {} atoms, {} "
							  "residues at {}>"),
				demangled_type<float>(), m_natoms, m_nresidues, fmt::ptr(this));
		}

		virtual float custom_kernel(const std::shared_ptr<Residue<float>>& residue,
			const std::shared_ptr<Residue<float>>& next_residue) const
		{
			throw "Error, custom kernel has not been defined";
		}

		arma::Col<float> run_custom_kernel() const
		{
			auto custom_lambda = [this](const std::shared_ptr<Residue<float>>& residue,
									 const std::shared_ptr<Residue<float>>& next_residue) -> float {
				return this->custom_kernel(residue, next_residue);
			};
			return arma::Col<float>(
				core::residue_kernel_engine(this->m_residues, 0, custom_lambda).memptr(),
				m_nresidues);
		}

		arma::Mat<float> run_custom_pairwise_kernel() const
		{
			auto custom_lambda = [this](const std::shared_ptr<Residue<float>>& residue,
									 const std::shared_ptr<Residue<float>>& next_residue) -> float {
				return this->custom_kernel(residue, next_residue);
			};
			return core::pairwise_residue_kernel_engine(this->m_residues, 0, custom_lambda).slice(0);
		}
	};
}
#endif // SWIGPYTHON

#endif // PROSTRUCT_CustomPDB_H
