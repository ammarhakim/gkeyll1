/**
 * @file	lctxyzfieldsetter.cxx
 *
 * @brief	Unit tests for TXYZFieldSetter updater
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcField.h>
#include <LcLuaTXYZFunction.h>
#include <LcRectCartGrid.h>
#include <LcTXYZFieldSetter.h>
#include <LcTest.h>

// std includes
#include <cmath>

void
test_1()
{
  int lower[2] = {0, 0};
  int upper[2] = {10, 12};
  int lg[2] = {0, 0};
  int ug[2] = {0, 0};
  Lucee::Region<2, int> rgn(lower, upper);
  Lucee::Field<2, double> a(rgn, 3, lg, ug, 10);

// create grid
  double lxy[2] = {0.0, 0.0};
  double uxy[2] = {1.0, 1.0};
  Lucee::Region<2, double> compSpace(lxy, uxy);
  Lucee::RectCartGrid<2> grid2(rgn, rgn, compSpace);

  Lucee::TXYZFieldSetter<2> txyzFieldSetter;
// call initialization sequence
  txyzFieldSetter.declareTypes();
  txyzFieldSetter.setGrid(grid2);
  txyzFieldSetter.initialize();

  std::vector<Lucee::DataStructIfc*> outVars(1);
// set input, output data-structures
  outVars[0] = &a;
  txyzFieldSetter.setOutVars(outVars);

// string with table
  std::string tblStr = 
    "func = {"
    "  f = function (t, x, y, z)"
    "    return t*x, t*y, t*z"
    "  end"
    "}";

  Lucee::LuaState L;
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "func");

// push name of function on stack
  lua_pushstring(L, "f");
  lua_gettable(L, -2);
  Lucee::LuaTXYZFunction fun(L, "f", 3);

  txyzFieldSetter.setFunObj(fun);
  txyzFieldSetter.setCurrTime(0.0);
// run updater
  txyzFieldSetter.update(1.0);

  double xc[3];
  Lucee::ConstFieldPtr<double> aPtr = a.createConstPtr();
  for (int i=a.getLower(0); i<a.getUpper(0); ++i)
    for (int j=a.getLower(1); j<a.getUpper(1); ++j)
    {
      grid2.setIndex(i, j);
      grid2.getCentriod(xc);
      a.setPtr(aPtr, i, j);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[0] == xc[0]);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[1] == xc[1]);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[2] == xc[2]);
    }

// run updater
  txyzFieldSetter.update(0.25);

  for (int i=a.getLower(0); i<a.getUpper(0); ++i)
    for (int j=a.getLower(1); j<a.getUpper(1); ++j)
    {
      grid2.setIndex(i, j);
      grid2.getCentriod(xc);
      a.setPtr(aPtr, i, j);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[0] == 0.25*xc[0]);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[1] == 0.25*xc[1]);
      LC_ASSERT("Testing txyzFieldSetter", aPtr[2] == 0.25*xc[2]);
    }
}

int
main(void)
{
  LC_BEGIN_TESTS("lctxyzfieldsetter");
  test_1();
  LC_END_TESTS;
}
