/**
 * @file	lcstructgridfield.cxx
 *
 * @brief	Unit tests for Lucee::StructGridFields class
 */

// lucee includes
#include <LcCartProdDecompRegionCalc.h>
#include <LcDecompRegion.h>
#include <LcGlobals.h>
#include <LcLuaState.h>
#include <LcRectCartGrid.h>
#include <LcRowMajorSequencer.h>
#include <LcStructGridField.h>
#include <LcTest.h>
#include <LcVec3.h>

// loki includes
#include <loki/Singleton.h>

// txbase includes
#ifdef HAVE_MPI
# include <TxMpiBase.h>
#else
# include <TxSelfBase.h>
#endif

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

// run some basic tests on grid
  LC_ASSERT("Testing grid", grid.getNumCells(0) == 10);
  LC_ASSERT("Testing grid", grid.getNumCells(1) == 15);

  LC_ASSERT("Testing grid spacing", grid.getDx(0) == 1.0/10.0);
  LC_ASSERT("Testing grid spacing", grid.getDx(1) == 1.5/15.0);

  double xc[3];
  int idx[2];
  Lucee::RowMajorSequencer<2> seq(grid.getLocalRegion());
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
// check cell volume
    LC_ASSERT("Checking cell volume", grid.getVolume() == grid.getDx(0)*grid.getDx(1));
// check cell centroid
    grid.getCentroid(xc);
    LC_ASSERT("Checking cell centroid", xc[0] == (idx[0]+0.5)*grid.getDx(0));
    LC_ASSERT("Checking cell centroid", xc[1] == (idx[1]+0.5)*grid.getDx(1));
    LC_ASSERT("Checking cell centroid", xc[2] == 0.0);
  }
}

void
test_1()
{
// create grid
  int lo[2] = {0, 0};
  int up[2] = {10, 15};
  Lucee::Region<2, int> localBox(lo, up);

  double plo[2] = {0.0, 0.0};
  double pup[2] = {1.0, 1.5};
  Lucee::Region<2, double> physBox(plo, pup);

  Lucee::RectCartGrid<2> grid(localBox, physBox);
// create field
  int lg[2] = {2, 3}, ug[2] = {3, 2};
  Lucee::StructGridField<2, double> fld(&grid, 3, lg, ug);
// check name
  LC_ASSERT("Testing name of field", fld.getName() == "Field2D");

// initialize field
  fld = 10.5;

// string with table
  std::string fnStr = 
    "  tbl = {init = function (x, y, z)"
    "    return x+10*y, x+20*y, x+30*y"
    "  end"
    " }";

  Lucee::LuaState L;
// evaluate string as Lua code
  if (luaL_loadstring(L, fnStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "tbl");

// push name of function on stack
  lua_pushstring(L, "init");
  lua_gettable(L, -2);
  if (! lua_isfunction(L, -1) )
  {
    lua_pop(L, 1);
    Lucee::Except lce("LuaTable::getFunctionRef: ");
    throw lce;
  }
  int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
  lua_pop(L, 1);
  fld.setFromLuaFunction(L, fnRef);

// test it
  double xc[3];
  int idx[2];
  Lucee::RowMajorSequencer<2> seq(grid.getLocalRegion());
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentroid(xc);
    LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 0) == xc[0]+10*xc[1]);
    LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 1) == xc[0]+20*xc[1]);
    LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 2) == xc[0]+30*xc[1]);
  }
}

void
test_2()
{
// create grid
  int lo[2] = {0, 0};
  int up[2] = {10, 15};
  Lucee::Region<2, int> localBox(lo, up);

  double plo[2] = {0.0, 0.0};
  double pup[2] = {1.0, 1.5};
  Lucee::Region<2, double> physBox(plo, pup);

  Lucee::RectCartGrid<2> grid(localBox, physBox);
// create field
  int lg[2] = {2, 3}, ug[2] = {3, 2};
  Lucee::StructGridField<2, double> fld(&grid, 3, lg, ug);
// check name
  LC_ASSERT("Testing name of field", fld.getName() == "Field2D");

// initialize field
  fld = 10.5;

// string with table
  std::string fnStr = 
    "  tbl = {init = function (x, y, z)"
    "    return x+10*y, x+20*y, x+30*y"
    "  end"
    " }";

  Lucee::LuaState L;
// evaluate string as Lua code
  if (luaL_loadstring(L, fnStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "tbl");

// push name of function on stack
  lua_pushstring(L, "init");
  lua_gettable(L, -2);
  if (! lua_isfunction(L, -1) )
  {
    lua_pop(L, 1);
    Lucee::Except lce("LuaTable::getFunctionRef: ");
    throw lce;
  }
  int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
  lua_pop(L, 1);
  fld.setGhostFromLuaFunction(L, fnRef, 0, 0); // x-direction, lower

// test it
  double xc[3];
  int idx[2];
  Lucee::RowMajorSequencer<2> seq(fld.getExtRegion());
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentroid(xc);
    if (idx[0] < fld.getLower(0))
    {
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 0) == xc[0]+10*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 1) == xc[0]+20*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 2) == xc[0]+30*xc[1]);
    }
    else
    {
      for (unsigned k=0; k<3; ++k)
        LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], k) == 10.5);
    }
  }

  fld = 10.5;
  fld.setGhostFromLuaFunction(L, fnRef, 0, 1); // x-direction, upper
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentroid(xc);
    if (idx[0] >= fld.getUpper(0))
    {
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 0) == xc[0]+10*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 1) == xc[0]+20*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 2) == xc[0]+30*xc[1]);
    }
    else
    {
      for (unsigned k=0; k<3; ++k)
        LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], k) == 10.5);
    }
  }

  fld = 10.5;
  fld.setGhostFromLuaFunction(L, fnRef, 1, 0); // y-direction, lower
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentroid(xc);
    if (idx[1] < fld.getLower(1))
    {
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 0) == xc[0]+10*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 1) == xc[0]+20*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 2) == xc[0]+30*xc[1]);
    }
    else
    {
      for (unsigned k=0; k<3; ++k)
        LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], k) == 10.5);
    }
  }

  fld = 10.5;
  fld.setGhostFromLuaFunction(L, fnRef, 1, 1); // y-direction, upper
  seq.reset();
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentroid(xc);
    if (idx[1] >= fld.getUpper(1))
    {
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 0) == xc[0]+10*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 1) == xc[0]+20*xc[1]);
      LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], 2) == xc[0]+30*xc[1]);
    }
    else
    {
      for (unsigned k=0; k<3; ++k)
        LC_ASSERT("Testing setting from Lua function", fld(idx[0], idx[1], k) == 10.5);
    }
  }
}

void
test_3()
{
// get communicator object
  TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
    ::Instance().comm;

// global region of grid
  int lo[2] = {0, 0};
  int up[2] = {64, 32};
  Lucee::Region<2, int> globalRgn(lo, up);

// create default decomp
  Lucee::DecompRegion<2> dcomp(globalRgn);

  int cuts[2] = {2, 2};
// create product decomposer
  Lucee::CartProdDecompRegionCalc<2> cartDecomp(cuts);

  if (comm->getNumProcs() == 4)
  {
    if (comm->getRank() == 0)
      std::cout << "Testing parallel field" << std::endl;
    cartDecomp.calcDecomp(comm->getNumProcs(), dcomp); // decompose

    double plo[2] = {0.0, 0.0};
    double pup[2] = {32, 32};
    Lucee::Region<2, double> physBox(plo, pup);
// create grid
    Lucee::RectCartGrid<2> grid(dcomp, physBox);

// create field
    int lg[2] = {2, 3}, ug[2] = {3, 2};
    Lucee::StructGridField<2, double> fld(&grid, 3, lg, ug);

// initialize field
    fld = 10.5;

// test it
    int idx[2];
    Lucee::ConstFieldPtr<double> cnstPtr = fld.createConstPtr();
    Lucee::RowMajorSequencer<2> seq(fld.getExtRegion());
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fld.setPtr(cnstPtr, idx);
      for (unsigned k=0; k<3; ++k)
        LC_ASSERT("Testing parallel field", cnstPtr[k] == 10.5);
    }

// now initialize interior of field
    Lucee::FieldPtr<double> ptr = fld.createPtr();
    Lucee::RowMajorSequencer<2> seqInt(fld.getRegion());
    while (seqInt.step())
    {
      seqInt.fillWithIndex(idx);
      fld.setPtr(ptr, idx);
      for (unsigned k=0; k<3; ++k)
        ptr[k] = 12.5;
    }

    Lucee::Region<2, int> rgn = fld.getRegion();
    seq.reset();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fld.setPtr(cnstPtr, idx);
      if (rgn.isInside(idx))
        for (unsigned k=0; k<3; ++k)
// ensure interior id 12.5
          LC_ASSERT("Testing parallel field interior", cnstPtr[k] == 12.5);
      else
        for (unsigned k=0; k<3; ++k)
// ensure ghosts are still 10.5
          LC_ASSERT("Testing parallel field ghost", cnstPtr[k] == 10.5);
    }

// now sync ghost cells values
    fld.sync();

    Lucee::Region<2, int> globalRegion = fld.getGlobalRegion();
// check if ghost values were updated correctly
    seq.reset();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      fld.setPtr(cnstPtr, idx);
      if (globalRegion.isInside(idx))
        for (unsigned k=0; k<3; ++k)
          LC_ASSERT("Testing parallel after sync", cnstPtr[k] == 12.5);
    }
  }
}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcstructgridfield");
  test_0();
  test_1();
  test_2();
  test_3();
  LC_END_TESTS;
}
