/**
 * @file	lcdeocmpregion.cxx
 *
 * @brief	Unit tests for Lucee::DecompRegion class
 */

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
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
}

void
test_2()
{
  int lower[2] = {0, 0};
  int upper[2] = {32, 64};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  unsigned cuts[2] = {2, 4};

// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

// first try to decompose with incorrect number of regions
  LC_RAISES("Testing if incorrect decomposition can be done", 
    cartDecomp.calcDecomp(10, dcomp), Lucee::Except);

// now decompose domain
  cartDecomp.calcDecomp(8, dcomp);
}

int
main(void) 
{
  LC_BEGIN_TESTS("lcdecomregion");
  test_1();
  test_2();
  LC_END_TESTS;
}
