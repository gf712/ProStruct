//
// Created by gil on 05/04/18.
//

#define BOOST_TEST_MODULE PDBTestSuite
#define BOOST_TEST_DYN_LINK

#include <boost/test/included/unit_test.hpp>
#include "PDB.h"

namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(PDBTests)

    BOOST_AUTO_TEST_CASE(LoadPDB) {

        PDB pdb = PDB("test.pdb");

        BOOST_TEST(pdb.n_chains() == 2);

        BOOST_TEST(pdb.getChainIDs()[0] == "L");
        BOOST_TEST(pdb.getChainIDs()[1] == "H");

        BOOST_TEST(pdb.getChain("H")->n_residues() == 123);
        BOOST_TEST(pdb.getChain("H")->n_atoms() == 974);

        BOOST_TEST(pdb.getChain("L")->n_residues() == 113);
        BOOST_TEST(pdb.getChain("L")->n_atoms() == 893);

        BOOST_TEST(pdb.getXYZ().n_rows == 3);
        BOOST_TEST(pdb.getXYZ().n_cols == 1867);

    }

    BOOST_AUTO_TEST_CASE(PredictBackboneHBonds) {

        PDB pdb = PDB("test.pdb");

        arma::mat E = pdb.predict_backboneHbonds();

        int total = std::accumulate(E.begin(), E.end(), 0);

        BOOST_TEST(E(2, 25) == 1);
        BOOST_TEST(total == 159);

    }

    BOOST_AUTO_TEST_CASE(ShrakeRupley) {

        PDB pdb = PDB("test.pdb");

        arma::vec asa = pdb.calculate_ASA(1.4);

        BOOST_TEST(asa.at(0) == 43.953897152668652, tt::tolerance(10e-9));

    }

    BOOST_AUTO_TEST_CASE(KabschSander) {

        PDB pdb = PDB("test.pdb");

        arma::mat E = pdb.calculate_KabschSander();

        BOOST_TEST(E.n_rows == pdb.n_residues());
        BOOST_TEST(E.n_cols == pdb.n_residues());

        BOOST_TEST(E(2, 25) == -1.9515432974890459, tt::tolerance(10e-9));

    }

    BOOST_AUTO_TEST_CASE(Kabsch_RMSD) {

        PDB pdb = PDB("test.pdb");

        double rmsd = pdb.kabsch_rmsd(pdb);

        BOOST_TEST(rmsd == 0.0, tt::tolerance(10e-7));
    }

BOOST_AUTO_TEST_SUITE_END()
