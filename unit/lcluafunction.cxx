/**
 * @file	lcluatable.cxx
 *
 * @brief	Unit tests for using LUA into Lucee
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcLuaFunction.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>
#include <LcTest.h>

void
test_1(Lucee::LuaState& L)
{
// string with table
  std::string tblStr = 
    "background = {"
    "  init = function (t, x, y, z)"
    "    return t*(x+y+z)"
    "  end"
    "}";

// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// push name of function on stack
  lua_pushstring(L, "init");
  lua_gettable(L, -2);
  Lucee::LuaFunction fun(L, "init", 1);

  std::vector<double> txyz(4);
  double t = 10.0;

  double dx=0.1, dy=0.2, dz=0.3;
  for (unsigned i=0; i<100; ++i)
    for (unsigned j=0; j<100; ++j)
      for (unsigned k=0; k<10; ++k)
      {
        double x = 1.0 + dx*i;
        double y = 2.0 + dy*j;
        double z = 3.0 + dz*k;
        txyz[0] = t; txyz[1] = x; txyz[2] = y; txyz[3] = z;
        std::vector<double> res = fun.eval(txyz);
        LC_ASSERT("Testing size of return vector", res.size() == 1);
        LC_ASSERT("Testing return results", res[0] == t*(x+y+z));
      }

  t = 10.0;
  for (unsigned i=0; i<100; ++i)
    for (unsigned j=0; j<100; ++j)
    {
      double x = 1.0 + dx*i;
      double y = 2.0 + dy*j;
      double z = 0.0;
      txyz[0] = t; txyz[1] = x; txyz[2] = y; txyz[3] = z;
      std::vector<double> res = fun.eval(txyz);
      LC_ASSERT("Testing size of return vector", res.size() == 1);
      LC_ASSERT("Testing return results", res[0] == t*(x+y+z));
    }
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluafunction");
  Lucee::LuaState L;
  test_1(L);

  LC_END_TESTS;
}
