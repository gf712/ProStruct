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

        std::string L = "L";
        std::string H = "H";

        BOOST_TEST(pdb.n_chains() == 2);

        BOOST_TEST(pdb.getChainIDs()[0] == L);
        BOOST_TEST(pdb.getChainIDs()[1] == H);

        BOOST_TEST(pdb.getChain("H")->n_residues() == 123);
        BOOST_TEST(pdb.getChain("H")->n_atoms() == 974);

        BOOST_TEST(pdb.getChain("L")->n_residues() == 113);
        BOOST_TEST(pdb.getChain("L")->n_atoms() == 893);

        BOOST_TEST(pdb.getXYZ().n_rows == 3);
        BOOST_TEST(pdb.getXYZ().n_cols == 1867);

    }

BOOST_AUTO_TEST_SUITE_END()
