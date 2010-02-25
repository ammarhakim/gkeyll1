/**
 * @file	lcrowmajorsequencer.cxx
 *
 * @brief	Unit tests for Lucee::RowMajorSequencer class
 *
 * @version	$Id: lcrowmajorsequencer.cxx 142 2009-08-23 16:32:29Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcRowMajorSequencer.h>
#include <LcTest.h>

void
test_1()
{
  int lower[1] = {3};
  int upper[1] = {12};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::RowMajorSequencer<1> seq(rgn);

  int idx[1];
  for (unsigned i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
  {
    seq.step();
    seq.fillWithIndex(idx);
    LC_ASSERT("Testing if 1D sequencer worked", idx[0] == i);
  }

  for (unsigned i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
  {
    seq.step();
    const int *cidx = seq.getIndex();
    LC_ASSERT("Testing if 1D sequencer worked", cidx[0] == i);
  }
}

void
test_2()
{
  int lower[2] = {2, -1};
  int upper[2] = {12, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::RowMajorSequencer<2> seq(rgn);

  int idx[2];
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    {
      seq.step();
      seq.fillWithIndex(idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j);
    }

  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    {
      seq.step();
      const int *cidx = seq.getIndex();
      LC_ASSERT("Testing if 2D sequencer worked", cidx[0] == i && cidx[1] == j);
    }
}

void
test_3()
{
  int lower[3] = {2, -1, 2};
  int upper[3] = {12, 10, 5};
  Lucee::Region<3, int> rgn(lower, upper);
  Lucee::RowMajorSequencer<3> seq(rgn);

  int idx[3];
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      for (int k=rgn.getLower(2); k<rgn.getUpper(2); ++k)
      {
        seq.step();
        seq.fillWithIndex(idx);
        LC_ASSERT("Testing if 3D sequencer worked", 
          idx[0] == i && idx[1] == j && idx[2] == k);
      }

  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
      for (int k=rgn.getLower(2); k<rgn.getUpper(2); ++k)
      {
        seq.step();
        const int *cidx = seq.getIndex();
        LC_ASSERT("Testing if 3D sequencer worked", 
          cidx[0] == i && cidx[1] == j && cidx[2] == k);
      }
}

void
test_4()
{
  int lower[2] = {2, -1};
  int upper[2] = {12, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::RowMajorSequencer<2> seq(rgn);

  int idx[2];
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    {
      seq.step();
      seq.fillWithIndex(idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j);
    }

  seq.reset();
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    {
      seq.step();
      seq.fillWithIndex(idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j);
    }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcrowmajorsequencer");
  test_1();
  test_2();
  test_3();
  test_4();
  LC_END_TESTS;
}
