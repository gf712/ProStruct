//
// Created by gil on 05/04/18.
//

#define BOOST_TEST_MODULE PDBTestSuite
#define BOOST_TEST_DYN_LINK

#include <boost/test/included/unit_test.hpp>
#include "PDB.h"

BOOST_AUTO_TEST_SUITE(PDBTests)

    BOOST_AUTO_TEST_CASE(LoadPDB) {

        auto pdb = PDB("test.pdb");

        BOOST_TEST(pdb.n_chains() == 2);

        BOOST_TEST(pdb.getChainIDs()[0] == "H");
        BOOST_TEST(pdb.getChainIDs()[1] == "L");

        BOOST_TEST(pdb.getChain("H")->n_residues() == 123);
        BOOST_TEST(pdb.getChain("H")->n_atoms() == 974);

        BOOST_TEST(pdb.getChain("L")->n_residues() == 113);
        BOOST_TEST(pdb.getChain("L")->n_atoms() == 893);

    }

BOOST_AUTO_TEST_SUITE_END()
