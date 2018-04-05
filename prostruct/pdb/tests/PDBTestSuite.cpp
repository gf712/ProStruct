//
// Created by gil on 05/04/18.
//

#define BOOST_TEST_MODULE PDBTestSuite
#define BOOST_TEST_DYN_LINK

#include <boost/test/included/unit_test.hpp>
#include "PDB.h"

BOOST_AUTO_TEST_SUITE(PDBTests)

    BOOST_AUTO_TEST_CASE(LoadPDB) {

        auto atom = PDB("test.pdb");

    }

BOOST_AUTO_TEST_SUITE_END()
