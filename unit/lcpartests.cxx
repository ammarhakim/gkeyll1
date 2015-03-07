/**
 * @file	lcpartests.cxx
 *
 * @brief	Testing various parallel ideas/things
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// txbase includes
#ifdef HAVE_MPI
# include <TxMpiBase.h>
#else
# include <TxSelfBase.h>
#endif

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcDecompRegion.h>
#include <LcGlobals.h>
#include <LcRectCartGrid.h>
#include <LcTest.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <iostream>
#include <vector>

bool isValid(TxCommBase *c)
{ return (bool) c; }

void
test_0()
{
// create grid
  int lo[2] = {0, 0};
  int up[2] = {10, 15};
  Lucee::Region<2, int> localBox(lo, up);

  double plo[2] = {0.0, 0.0};
  double pup[2] = {1.0, 1.5};
  Lucee::Region<2, double> physBox(plo, pup);

  Lucee::RectCartGrid<2> grid(localBox, physBox);

  std::vector<unsigned> cdirs(1);
  cdirs[0] = 0;
// create sub-grid
  //Lucee::RectCartGrid<1> subGrid = grid.createSubGrid<1>(cdirs);
// check it
  //LC_ASSERT("Testing grid", subGrid.getNumCells(0) == 10);
  //LC_ASSERT("Testing grid spacing", subGrid.getDx(0) == 1.0/10.0);
}

int
main(int argc, char *argv[])
{
  LC_BEGIN_TESTS("lcpartests");

// get hold of global communicator
  TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
    ::Instance().comm;

// create a decomposition
  int lower[2] = {0, 0};
  int upper[2] = {64, 128};
  Lucee::Region<2, int> globalRgn(lower, upper);
// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);
  int cuts[2] = {2, 2};
// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

  int numProcs = comm->getNumProcs();
  cartDecomp.calcDecomp(numProcs, dcomp);

  test_0();

  LC_END_TESTS;
  return 0;
}
