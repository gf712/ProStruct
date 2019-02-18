//
// Created by gil on 27/03/18.
//

#include "gtest/gtest.h"

#include "prostruct/struct/bond.h"

TEST(BondTest, CreateBond1)
{
	auto bond = Bond<double>(1., 1., 1., 1);

	ASSERT_EQ(bond.getBondVector()[0], 1.);
	ASSERT_EQ(bond.getBondVector()[1], 1.);
	ASSERT_EQ(bond.getBondVector()[2], 1.);
	ASSERT_EQ(bond.getBondType(), 1);
}