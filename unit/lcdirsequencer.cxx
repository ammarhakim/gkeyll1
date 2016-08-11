/**
 * @file	lcrowmajorsequencer.cxx
 *
 * @brief	Unit tests for Lucee::RowMajorSequencer class
 */

// lucee includes
#include <LcDirSequencer.h>
#include <LcTest.h>

void
test_1()
{
  int lower[1] = {3};
  int upper[1] = {12};
  Lucee::Region<1, int> rgn(lower, upper);
  Lucee::DirSequencer<1> seq(rgn, 0, 1, 2);

  int idx[1];
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
  {
    seq.step();
    seq.fillWithIndex(0, idx);
    LC_ASSERT("Testing if 1D sequencer worked", idx[0] == i);

    seq.fillWithIndex(-1, idx); // i-1
    LC_ASSERT("Testing if 1D sequencer worked", idx[0] == i-1);

    seq.fillWithIndex(1, idx); // i+1
    LC_ASSERT("Testing if 1D sequencer worked", idx[0] == i+1);

    seq.fillWithIndex(2, idx); // i+2
    LC_ASSERT("Testing if 1D sequencer worked", idx[0] == i+2);
  }
}

void
test_2()
{
  int lower[2] = {2, -1};
  int upper[2] = {12, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::DirSequencer<2> seq(rgn, 0, 1, 2);

  int idx[2];
  for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    {
      seq.step();
      seq.fillWithIndex(0, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j);

      seq.fillWithIndex(-1, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i-1 && idx[1] == j);

      seq.fillWithIndex(1, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i+1 && idx[1] == j);

      seq.fillWithIndex(2, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i+2 && idx[1] == j);
    }
}

void
test_3()
{
  int lower[2] = {2, -1};
  int upper[2] = {12, 10};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::DirSequencer<2> seq(rgn, 1, 1, 2);

  int idx[2];
  for (int i=rgn.getLower(0); i<rgn.getUpper(0); ++i)
    for (int j=rgn.getLower(1); j<rgn.getUpper(1); ++j)
    {
      seq.step();
      seq.fillWithIndex(0, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j);

      seq.fillWithIndex(-1, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j-1);

      seq.fillWithIndex(1, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j+1);

      seq.fillWithIndex(2, idx);
      LC_ASSERT("Testing if 2D sequencer worked", idx[0] == i && idx[1] == j+2);
    }
}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcdirsequencer");
  test_1();
  test_2();
  test_3();
  LC_END_TESTS;
}
