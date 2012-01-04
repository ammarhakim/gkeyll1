/**
 * @file	lcdeocmpregion.cxx
 *
 * @brief	Unit tests for Lucee::DecompRegion class
 */

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcColMajorSequencer.h>
#include <LcDecompRegion.h>
#include <LcTest.h>

void
test_1()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);

// run tests on default decomp (which is unitary)
  LC_ASSERT("Testing unitary decomp", dcomp.getNumRegions() == 1);
  Lucee::Region<2, int> subRgn = dcomp.getRegion(0);
  for (unsigned i=0; i<2; ++i)
  {
    LC_ASSERT("Testing unitary decomp", subRgn.getLower(i) == globalRgn.getLower(i));
    LC_ASSERT("Testing unitary decomp", subRgn.getUpper(i) == globalRgn.getUpper(i));
  }
  LC_ASSERT("Testing unitary decomp", dcomp.checkCovering() == true);

  LC_ASSERT("Tesing unitary decomp", dcomp.calcMinMaxVolRatio() == 1.0);
}

void
test_2()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {2, 4};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// first try to decompose with incorrect number of regions
  LC_RAISES("Testing if incorrect decomposition can be done", 
    cartDecomp.calcDecomp(10, dcomp), Lucee::Except);

// now decompose domain
  cartDecomp.calcDecomp(8, dcomp);
// check decomp
  int xlower[2] = {0, 16};
  int ylower[4] = {0, 16, 32, 48};

  LC_ASSERT("Testing decomposition", dcomp.getNumRegions() == 8);

  Lucee::Region<2, int> cutRgn(cuts);
  Lucee::ColMajorSequencer<2> seq(cutRgn);
  unsigned r = 0;
  int idx[2];
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    Lucee::Region<2, int> subRgn = dcomp.getRegion(r);
    LC_ASSERT("Testing decomposition", subRgn.getLower(0) == xlower[idx[0]]);
    LC_ASSERT("Testing decomposition", subRgn.getLower(1) == ylower[idx[1]]);

    LC_ASSERT("Testing decomposition", subRgn.getShape(0) == 16);
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == 16);

    r++;
  }

  LC_ASSERT("Tesing decomposition", dcomp.calcMinMaxVolRatio() == 1.0);
}

void
test_3()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {2, 8};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// now decompose domain
  cartDecomp.calcDecomp(16, dcomp);
// check decomp
  int xlower[2] = {0, 16};
  int ylower[8] = {0, 8, 16, 24, 32, 40, 48, 56};

  LC_ASSERT("Testing decomposition", dcomp.getNumRegions() == 16);

  Lucee::Region<2, int> cutRgn(cuts);
  Lucee::ColMajorSequencer<2> seq(cutRgn);
  unsigned r = 0;
  int idx[2];
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    Lucee::Region<2, int> subRgn = dcomp.getRegion(r);
    LC_ASSERT("Testing decomposition", subRgn.getLower(0) == xlower[idx[0]]);
    LC_ASSERT("Testing decomposition", subRgn.getLower(1) == ylower[idx[1]]);

    LC_ASSERT("Testing decomposition", subRgn.getShape(0) == 16);
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == 8);

    r++;
  }

  LC_ASSERT("Tesing decomposition", dcomp.calcMinMaxVolRatio() == 1.0);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcdecomregion");
  test_1();
  test_2();
  test_3();
  LC_END_TESTS;
}
