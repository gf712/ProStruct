//
// Created by gil on 05/04/18.
//

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

TYPED_TEST(PDBTest, LoadPDB) {

    auto pdb = PDB<TypeParam>("test.pdb");

    ASSERT_EQ(pdb.n_chains(), 2);

    ASSERT_EQ(pdb.getChainIDs()[0], "L");
    ASSERT_EQ(pdb.getChainIDs()[1], "H");

    ASSERT_EQ(pdb.getChain("H")->n_residues(), 123);
    ASSERT_EQ(pdb.getChain("H")->n_atoms(), 974);

    ASSERT_EQ(pdb.getChain("L")->n_residues(), 113);
    ASSERT_EQ(pdb.getChain("L")->n_atoms(), 893);

    ASSERT_EQ(pdb.getXYZ().n_rows, 3);
    ASSERT_EQ(pdb.getXYZ().n_cols, 1867);

}

TYPED_TEST(PDBTest, PredictBackboneHBonds) {

    auto pdb = PDB<TypeParam>("test.pdb");

    auto E = pdb.predict_backboneHbonds();

    int total = std::accumulate(E.begin(), E.end(), 0);

    ASSERT_EQ(E(2, 25), 1);
    ASSERT_EQ(total, 159);

}

TYPED_TEST(PDBTest, ShrakeRupley) {

    auto pdb = PDB<TypeParam>("test.pdb");
    auto asa = pdb.calculate_ASA(1.4);

    EXPECT_NEAR(asa.at(0), 43.953897152668652, get_epsilon<TypeParam>());
}

TYPED_TEST(PDBTest, KabschSander) {

    auto pdb = PDB<TypeParam>("test.pdb");

    auto E = pdb.calculate_KabschSander();

    ASSERT_EQ(E.n_rows, pdb.n_residues());
    ASSERT_EQ(E.n_cols, pdb.n_residues());

    EXPECT_NEAR(E(2, 25), -1.9515432974890459, get_epsilon<TypeParam>());

}

TYPED_TEST(PDBTest, Kabsch_RMSD) {

    auto pdb = PDB<TypeParam>("test.pdb");

    auto rmsd = pdb.kabsch_rmsd(pdb);

    EXPECT_NEAR(rmsd, 0.0, get_epsilon<TypeParam>());
}
