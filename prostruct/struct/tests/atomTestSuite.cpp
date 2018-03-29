//
// Created by gil on 27/03/18.
//

#define BOOST_TEST_MODULE AtomTestSuite
#define BOOST_TEST_DYN_LINK

#include <boost/test/included/unit_test.hpp>
#include "atom.h"

BOOST_AUTO_TEST_SUITE(AtomTests)

    BOOST_AUTO_TEST_CASE(LoadAtom) {

        auto atom = Atom("H");

        BOOST_TEST(atom.getAtomicNumber(), 1);
        BOOST_TEST(atom.getAtomicWeight(), 1.0079);
    }

BOOST_AUTO_TEST_SUITE_END()
