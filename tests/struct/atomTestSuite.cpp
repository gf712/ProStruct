//
// Created by gil on 27/03/18.
//

#include "gtest/gtest.h"

#include "prostruct/struct/atom.h"

TEST(AtomTests, LoadAtom)
{
  auto atom = Atom("H");

  ASSERT_EQ(atom.getAtomicNumber(), 1);
  ASSERT_EQ(atom.getAtomicWeight(), 1.0079);
}