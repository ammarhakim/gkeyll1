/**
 * @file	lcstructgridfield.cxx
 *
 * @brief	Unit tests for Lucee::StructGridFields class
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLuaState.h>
#include <LcRectCartGrid.h>
#include <LcStructGridField.h>
#include <LcTest.h>

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

  Lucee::RectCartGrid<2> grid(localBox, localBox, physBox);

// run some basic tests on grid
  LC_ASSERT("Testing grid", grid.getNumCells(0) == 10);
  LC_ASSERT("Testing grid", grid.getNumCells(1) == 15);

  LC_ASSERT("Testing grid spacing", grid.getDx(0) == 1.0/10.0);
  LC_ASSERT("Testing grid spacing", grid.getDx(1) == 1.5/15.0);

  double xc[3];
  int idx[2];
  Lucee::RowMajorSequencer<2> seq(grid.getLocalBox());
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
// check cell volume
    LC_ASSERT("Checking cell volume", grid.getVolume() == grid.getDx(0)*grid.getDx(1));
// check cell centroid
    grid.getCentriod(xc);
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

  Lucee::RectCartGrid<2> grid(localBox, localBox, physBox);
// create field
  int lg[2] = {2, 3}, ug[2] = {3, 2};
  Lucee::StructGridField<2, double> fld(&grid, 3, lg, ug);
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
  Lucee::RowMajorSequencer<2> seq(grid.getLocalBox());
  while (seq.step())
  {
    seq.fillWithIndex(idx);
    grid.setIndex(idx);
    grid.getCentriod(xc);
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

  Lucee::RectCartGrid<2> grid(localBox, localBox, physBox);
// create field
  int lg[2] = {2, 3}, ug[2] = {3, 2};
  Lucee::StructGridField<2, double> fld(&grid, 3, lg, ug);
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
    grid.getCentriod(xc);
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
    grid.getCentriod(xc);
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
    grid.getCentriod(xc);
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
    grid.getCentriod(xc);
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

int
main(void)
{
  LC_BEGIN_TESTS("lcstructgridfield");
  test_0();
  test_1();
  test_2();
  LC_END_TESTS;
}
