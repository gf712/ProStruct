/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "prostruct/struct/chain.h"
#include <gtest/gtest.h>

using namespace prostruct;

TEST(ResidueTest, SimpleTwoAASequence)
{

	auto N1 = std::make_shared<Atom<double>>("N", "N", 35.446, 51.519, 6.329);
	auto CA1 = std::make_shared<Atom<double>>("C", "CA", 35.098, 51.281, 4.896);
	auto C1 = std::make_shared<Atom<double>>("C", "C", 33.588, 51.392, 4.684);
	auto O1 = std::make_shared<Atom<double>>("O", "O", 33.012, 50.680, 3.858);
	auto CB1 = std::make_shared<Atom<double>>("C", "CB", 35.595, 49.896, 4.462);
	auto CG1 = std::make_shared<Atom<double>>("C", "CG", 34.919, 48.765, 5.220);
	auto OD11 = std::make_shared<Atom<double>>("O", "OD1", 34.684, 48.919, 6.438);
	auto OD21 = std::make_shared<Atom<double>>("O", "OD2", 34.635, 47.715, 4.600);

	auto N2 = std::make_shared<Atom<double>>("N", "N", 32.964, 52.298, 5.433);
	auto CA2 = std::make_shared<Atom<double>>("C", "CA", 31.521, 52.533, 5.366);
	auto C2 = std::make_shared<Atom<double>>("C", "C", 31.011, 52.774, 3.947);
	auto O2 = std::make_shared<Atom<double>>("O", "O", 30.021, 52.173, 3.525);
	auto CB2 = std::make_shared<Atom<double>>("C", "CB", 31.143, 53.743, 6.231);
	auto CG2 = std::make_shared<Atom<double>>("C", "CG", 31.514, 53.618, 7.696);
	auto CD2 = std::make_shared<Atom<double>>("C", "CD", 31.171, 54.883, 8.482);
	auto NE2 = std::make_shared<Atom<double>>("N", "NE", 29.739, 55.187, 8.480);
	auto CZ2 = std::make_shared<Atom<double>>("C", "CZ", 29.141, 56.009, 7.622);
	auto NH12 = std::make_shared<Atom<double>>("N", "NH1", 29.847, 56.622, 6.682);
	auto NH22 = std::make_shared<Atom<double>>("N", "NH2", 27.832, 56.219, 7.704);

	auto asp = std::make_shared<Residue<double>>(atomVector<double>({ N1, CA1, C1, O1, CB1, CG1, OD11, OD21 }),
		"ASP", "ASP1");
	auto arg = std::make_shared<Residue<double>>(atomVector<double>({ N2, CA2, C2, O2, CB2, CG2, CD2, CD2,
													 NE2, CZ2, NH12, NH22 }),
		"ARG", "ARG2");

	auto chain = Chain<double>(residueVector<double>({ asp, arg }), "Chain1");

	ASSERT_EQ(asp->getBackbone()[2]->getName(), "C");
	ASSERT_EQ(arg->getBackbone()[0]->getName(), "N");

	ASSERT_TRUE(asp->getBackbone()[2]->hasBond(arg->getBackbone()[0]));
}