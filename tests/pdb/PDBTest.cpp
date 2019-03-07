/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "gtest/gtest.h"

#include "prostruct/prostruct.h"

using namespace prostruct;

template <typename T>
constexpr T get_epsilon()
{
	return std::numeric_limits<T>::epsilon();
}
template <>
constexpr float get_epsilon()
{
	return 0.01; // rmsd gets quite imprecise when working with floats
}
template <>
constexpr double get_epsilon()
{
	return 1e-05;
}

template <typename T>
class PDBTest : public ::testing::Test {
};
using floatTypes = ::testing::Types<float, double>;

TYPED_TEST_CASE(PDBTest, floatTypes);

TYPED_TEST(PDBTest, LoadPDB)
{

	auto pdb = PDB<TypeParam>("test.pdb");

	ASSERT_EQ(pdb.n_chains(), 2);

	ASSERT_EQ(pdb.get_chain_names()[0], "L");
	ASSERT_EQ(pdb.get_chain_names()[1], "H");

	ASSERT_EQ(pdb.get_chain("H")->n_residues(), 123);
	ASSERT_EQ(pdb.get_chain("H")->n_atoms(), 974);

	ASSERT_EQ(pdb.get_chain("L")->n_residues(), 113);
	ASSERT_EQ(pdb.get_chain("L")->n_atoms(), 893);

	ASSERT_EQ(pdb.get_xyz().n_rows, 3);
	ASSERT_EQ(pdb.get_xyz().n_cols, 1867);
}

TYPED_TEST(PDBTest, PredictBackboneHBonds)
{

	auto pdb = PDB<TypeParam>("test.pdb");

	auto E = pdb.predict_backbone_hbonds();

	int total = std::accumulate(E.begin(), E.end(), 0);

	ASSERT_EQ(E(2, 25), 1);
	ASSERT_EQ(total, 160);
}

TYPED_TEST(PDBTest, ShrakeRupley)
{

	auto pdb = PDB<TypeParam>("test.pdb");
	auto asa = pdb.compute_shrake_rupley(1.4, 960);

	EXPECT_NEAR(asa.at(0), 43.593459609528409, get_epsilon<TypeParam>());
}

TYPED_TEST(PDBTest, KabschSander)
{

	auto pdb = PDB<TypeParam>("test.pdb");

	auto E = pdb.compute_kabsch_sander();

	ASSERT_EQ(E.n_rows, pdb.n_residues());
	ASSERT_EQ(E.n_cols, pdb.n_residues());

	EXPECT_NEAR(E(2, 25), -1.9521032812185981, get_epsilon<TypeParam>());
}

TYPED_TEST(PDBTest, Kabsch_RMSD)
{

	auto pdb = PDB<TypeParam>("test.pdb");

	auto rmsd = pdb.kabsch_rmsd(pdb);

	EXPECT_NEAR(rmsd, 0.0, get_epsilon<TypeParam>());
}

TYPED_TEST(PDBTest, phi_angles)
{
	auto pdb = PDB<TypeParam>("test.pdb");

	auto phi_degree = pdb.calculate_phi();

	EXPECT_NEAR(phi_degree(10), -76.065502471, get_epsilon<TypeParam>());

	auto phi_rad = pdb.calculate_phi(true);

	EXPECT_NEAR(phi_rad(10), -1.3275937, get_epsilon<TypeParam>());
}

TYPED_TEST(PDBTest, psi_angles)
{
	auto pdb = PDB<TypeParam>("test.pdb");

	auto psi_degree = pdb.calculate_psi();

	EXPECT_NEAR(psi_degree(10), 94.52472, get_epsilon<TypeParam>());

	auto psi_rad = pdb.calculate_psi(true);

	EXPECT_NEAR(psi_rad(10), 1.6497685, get_epsilon<TypeParam>());
}
