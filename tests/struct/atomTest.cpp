/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE', which is part of this source code package.
 *
 * Authors: Gil Hoben
 *
 */

#include "gtest/gtest.h"

#include "prostruct/struct/atom.h"

TEST(AtomTest, LoadAtom)
{
	auto atom = Atom<double>("H");

	ASSERT_EQ(atom.getAtomicNumber(), 1);
	ASSERT_EQ(atom.getAtomicWeight(), 1.0079);
}