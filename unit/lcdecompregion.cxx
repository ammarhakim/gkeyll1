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
test_2_2()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {4, 2};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// first try to decompose with incorrect number of regions
  LC_RAISES("Testing if incorrect decomposition can be done", 
    cartDecomp.calcDecomp(10, dcomp), Lucee::Except);

// now decompose domain
  cartDecomp.calcDecomp(8, dcomp);
// check decomp
  int xlower[4] = {0, 8, 16, 24};
  int ylower[2] = {0, 32};

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

    LC_ASSERT("Testing decomposition", subRgn.getShape(0) == 8);
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == 32);

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

void
test_4()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 70};
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
  int ylower[8] = {0, 8+1, 16+2, 24+3, 32+4, 40+5, 48+6, 56+6};
  int shape[8] = {9, 9, 9, 9, 9, 9, 8, 8};

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
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == shape[idx[1]]);

    r++;
  }

  LC_ASSERT("Tesing decomposition", dcomp.calcMinMaxVolRatio() == 8.0/9.0);
}

void
test_5()
{
  int lower[2] = {0, 0};
  int upper[2] = {35, 70};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {2, 8};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// now decompose domain
  cartDecomp.calcDecomp(16, dcomp);
// check decomp
  int xlower[2] = {0, 18};
  int xshape[2] = {18, 17};

  int ylower[8] = {0, 8+1, 16+2, 24+3, 32+4, 40+5, 48+6, 56+6};
  int yshape[8] = {9, 9, 9, 9, 9, 9, 8, 8};

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

    LC_ASSERT("Testing decomposition", subRgn.getShape(0) == xshape[idx[0]]);
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == yshape[idx[1]]);

    r++;
  }

  LC_ASSERT("Tesing decomposition", dcomp.calcMinMaxVolRatio() == (8.0*17.0)/(9.0*18));
}

void
test_6()
{
  int lower[4] = {0, 0, 0, 0};
  int upper[4] = {32, 64, 32, 64};
  Lucee::Region<4, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<4> dcomp(globalRgn);
  int cuts[4] = {2, 8, 2, 8};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<4> cartDecomp(cuts);

// now decompose domain
  cartDecomp.calcDecomp(2*8*2*8, dcomp);
// check decomp
  int x0lower[2] = {0, 16};
  int x1lower[8] = {0, 8, 16, 24, 32, 40, 48, 56};
  int x2lower[2] = {0, 16};
  int x3lower[8] = {0, 8, 16, 24, 32, 40, 48, 56};

  LC_ASSERT("Testing decomposition", dcomp.getNumRegions() == 2*8*2*8);

  Lucee::Region<4, int> cutRgn(cuts);
  Lucee::ColMajorSequencer<4> seq(cutRgn);
  unsigned r = 0;
  int idx[4];
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    Lucee::Region<4, int> subRgn = dcomp.getRegion(r);
    LC_ASSERT("Testing decomposition", subRgn.getLower(0) == x0lower[idx[0]]);
    LC_ASSERT("Testing decomposition", subRgn.getLower(1) == x1lower[idx[1]]);
    LC_ASSERT("Testing decomposition", subRgn.getLower(2) == x2lower[idx[2]]);
    LC_ASSERT("Testing decomposition", subRgn.getLower(3) == x3lower[idx[3]]);

    LC_ASSERT("Testing decomposition", subRgn.getShape(0) == 16);
    LC_ASSERT("Testing decomposition", subRgn.getShape(1) == 8);
    LC_ASSERT("Testing decomposition", subRgn.getShape(2) == 16);
    LC_ASSERT("Testing decomposition", subRgn.getShape(3) == 8);

    r++;
  }

  LC_ASSERT("Tesing decomposition", dcomp.calcMinMaxVolRatio() == 1.0);
}

void
test_7()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {2, 4};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// now decompose domain
  cartDecomp.calcDecomp(8, dcomp);

  unsigned numNeigh[8] = {3, 3, 5, 5, 5, 5, 3, 3};

  int lowerExt[2] = {1, 1};
  int upperExt[2] = {1, 1};
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getNeighbors(rn, lowerExt, upperExt);

    LC_ASSERT("Testing number of neighbors", neigh.size() == numNeigh[rn]);
    Lucee::Region<2, int> rgn = dcomp.getRegion(rn);
// now check if these are really neighbors
    std::vector<unsigned>::const_iterator itr = neigh.begin();
    for ( ; itr!=neigh.end(); ++itr)
    {
      Lucee::Region<2, int> neighRgn = dcomp.getRegion(*itr);
      LC_ASSERT("Testing if neighbors are correct",
        rgn.extend(lowerExt, upperExt).intersect(neighRgn).isEmpty() == false);
    }
  }
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcdecomregion");
  test_1();
  test_2();
  test_2_2();
  test_3();
  test_4();
  test_5();
  test_6();
  test_7();
  LC_END_TESTS;
}
