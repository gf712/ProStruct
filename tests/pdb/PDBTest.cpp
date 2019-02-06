//
// Created by gil on 05/04/18.
//

#include "gtest/gtest.h"

#include "prostruct/prostruct.h"

TEST(PDBTest, LoadPDB) {

    auto pdb = PDB<double>("test.pdb");

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

TEST(PDBTest, PredictBackboneHBonds) {

    auto pdb = PDB<double>("test.pdb");

    auto E = pdb.predict_backboneHbonds();

    int total = std::accumulate(E.begin(), E.end(), 0);

    ASSERT_EQ(E(2, 25), 1);
    ASSERT_EQ(total, 159);

}

//TEST(PDBTest, ShrakeRupley) {
//
//    PDB pdb = PDB("test.pdb");
//
//    arma::Col<T> asa = pdb.calculate_ASA(1.4);
//
//    EXPECT_NEAR(asa.at(0), 43.953897152668652, 10e-9);
//
//}

TEST(PDBTest, KabschSander) {

    auto pdb = PDB<double>("test.pdb");

    auto E = pdb.calculate_KabschSander();

    ASSERT_EQ(E.n_rows, pdb.n_residues());
    ASSERT_EQ(E.n_cols, pdb.n_residues());

    EXPECT_NEAR(E(2, 25), -1.9515432974890459, 10e-9);

}

TEST(PDBTest, Kabsch_RMSD) {

    auto pdb = PDB<double>("test.pdb");

    auto rmsd = pdb.kabsch_rmsd(pdb);

    EXPECT_NEAR(rmsd, 0.0, 10e-7);
}
