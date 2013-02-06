/**
 * @file	lcdeocmpregion.cxx
 *
 * @brief	Unit tests for Lucee::DecompRegion class
 */

// lucee includes
#include <LcCartGeneralDecompRegionCalc.h>
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
    std::vector<unsigned> neigh = dcomp.getRecvNeighbors(rn, lowerExt, upperExt);

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

// create copy and check
  Lucee::DecompRegion<2> dcompCpy(dcomp);
  LC_ASSERT("Checking if decompositions are the same", dcomp.compareDecomp(dcompCpy));

  Lucee::DecompRegion<2> dcompCpyFake(globalRgn);
  LC_ASSERT("Checking if decompositions are the same", dcomp.compareDecomp(dcompCpyFake) == false);
}

void
test_8()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);

  int lowerExt[2] = {1, 1};
  int upperExt[2] = {1, 1};
// check default decomp neighbors
  std::vector<unsigned> neigh = dcomp.getRecvNeighbors(0, lowerExt, upperExt);

  LC_ASSERT("Testing default decomp neighbors", neigh.size() == 0);
}

void
test_9()
{
  int lower[2] = {0, 0};
  int upper[2] = {50, 50};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);

// create product decomposer
  Lucee::CartGeneralDecompRegionCalc<2> cartDecomp;
// now decompose domain
  cartDecomp.calcDecomp(4, dcomp);

// check if decomp is covering
  LC_ASSERT("Testing if general decomp is covering", dcomp.checkCovering() == true);
}

void
test_10()
{
  int lower[3] = {0, 0, 0};
  int upper[3] = {50, 50, 201};
  Lucee::Region<3, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<3> dcomp(globalRgn);

// create product decomposer
  Lucee::CartGeneralDecompRegionCalc<3> cartDecomp;
// now decompose domain
  cartDecomp.calcDecomp(7, dcomp);

// check if decomp is covering
  LC_ASSERT("Testing if general decomp is covering", dcomp.checkCovering() == true);
}

void
test_11()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 10};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {1, 1};
  std::vector<int> periodicDirs(2);
  periodicDirs[0] = 0;

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);
// set periodic directions
  cartDecomp.setPeriodicDir(0);

// now decompose domain
  cartDecomp.calcDecomp(1, dcomp);

  LC_ASSERT("Testing decomposition", dcomp.getNumRegions() == 1);
  LC_ASSERT("Testing decomposition", dcomp.getNumTotalRegions() == 3);

// check decomposition
  Lucee::Region<2, int> rgn = dcomp.getRegion(1); // left pseudo box
  LC_ASSERT("Checking lower", rgn.getLower(0) == -10);
  LC_ASSERT("Checking lower", rgn.getLower(1) == 0);

  LC_ASSERT("Checking upper", rgn.getUpper(0) == 0);
  LC_ASSERT("Checking upper", rgn.getUpper(1) == 10);

  rgn = dcomp.getRegion(2); // right pseudo box
  LC_ASSERT("Checking lower", rgn.getLower(0) == 10);
  LC_ASSERT("Checking lower", rgn.getLower(1) == 0);

  LC_ASSERT("Checking upper", rgn.getUpper(0) == 20);
  LC_ASSERT("Checking upper", rgn.getUpper(1) == 10);

  int lowerExt[2], upperExt[2];

  lowerExt[0] = 1; lowerExt[1] = 0;
  upperExt[0] = 1; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getRecvNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 2);
  }

  lowerExt[0] = 1; lowerExt[1] = 0;
  upperExt[0] = 0; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getRecvNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 1);
    LC_ASSERT("Testing periodic BCs", neigh[0] == 1);
  }

  lowerExt[0] = 0; lowerExt[1] = 0;
  upperExt[0] = 1; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getRecvNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 1);
    LC_ASSERT("Testing periodic BCs", neigh[0] == 2);
  }

  lowerExt[0] = 1; lowerExt[1] = 0;
  upperExt[0] = 1; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getSendNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 2);
  }

  lowerExt[0] = 1; lowerExt[1] = 0;
  upperExt[0] = 0; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getSendNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 1);
    LC_ASSERT("Testing periodic BCs", neigh[0] == 2);
  }

  lowerExt[0] = 0; lowerExt[1] = 0;
  upperExt[0] = 1; upperExt[1] = 0;
// loop over each region, getting neighbors
  for (unsigned rn=0; rn<dcomp.getNumRegions(); ++rn)
  {
    std::vector<unsigned> neigh = dcomp.getSendNeighbors(rn, lowerExt, upperExt);
    LC_ASSERT("Testing periodic BCs", neigh.size() == 1);
    LC_ASSERT("Testing periodic BCs", neigh[0] == 1);
  }
}

void
test_12()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 10};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {1, 1};
  std::vector<int> periodicDirs(2);
  periodicDirs[0] = 0;

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);
// set periodic directions
  cartDecomp.setPeriodicDir(1);

// now decompose domain
  cartDecomp.calcDecomp(1, dcomp);

  LC_ASSERT("Testing decomposition", dcomp.getNumRegions() == 1);
  LC_ASSERT("Testing decomposition", dcomp.getNumTotalRegions() == 3);

// check decomposition
  Lucee::Region<2, int> rgn = dcomp.getRegion(1); // bottom pseudo box
  LC_ASSERT("Checking lower", rgn.getLower(0) == 0);
  LC_ASSERT("Checking lower", rgn.getLower(1) == -10);

  LC_ASSERT("Checking upper", rgn.getUpper(0) == 10);
  LC_ASSERT("Checking upper", rgn.getUpper(1) == 0);

  rgn = dcomp.getRegion(2); // top pseudo box
  LC_ASSERT("Checking lower", rgn.getLower(0) == 0);
  LC_ASSERT("Checking lower", rgn.getLower(1) == 10);

  LC_ASSERT("Checking upper", rgn.getUpper(0) == 10);
  LC_ASSERT("Checking upper", rgn.getUpper(1) == 20);
}

int
main(int argc, char **argv) 
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
  test_8();
  test_9();
  test_10();
  test_11();
  test_12();
  LC_END_TESTS;
}
