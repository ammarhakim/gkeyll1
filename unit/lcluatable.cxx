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
    "__kind = \"color\","
    "__type = \"Lucee\","
    "r = 0.3, b = 0.1, g = 0.0,"
    "name = \"blue_green\","
    "cells = {100, 50},"
    "address = {\"hello\", \"world\"},"
    "subcolors = {a=1.0, b=2.0, c=3.0},"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// construct LuaTable object
  Lucee::LuaTable back(L, "background");

  LC_ASSERT("Testing kind field", back.getKind() == "color");
  LC_ASSERT("Testing type field", back.getType() == "Lucee");

// test it
  LC_ASSERT("Testing Lua table", back.hasNumber("r"));
  LC_ASSERT("Testing Lua table", back.hasNumber("b"));
  LC_ASSERT("Testing Lua table", back.hasNumber("g"));
  LC_ASSERT("Testing Lua table", back.hasString("name"));
  LC_ASSERT("Testing Lua table", back.hasNumber("zzz") == false);

  LC_ASSERT("Testing Lua table", back.getNumber("r") == 0.3);
  LC_ASSERT("Testing Lua table", back.getNumber("b") == 0.1);
  LC_ASSERT("Testing Lua table", back.getNumber("g") == 0.0);
  LC_ASSERT("Testing Lua table", back.getString("name") == "blue_green");

// get subtable
  LC_ASSERT("Testing Lua table", back.hasTable("subcolors"));
  Lucee::LuaTable sub = back.getTable("subcolors");

  LC_ASSERT("Testing Lua table", sub.getNumber("a") == 1.0);
  LC_ASSERT("Testing Lua table", sub.getNumber("b") == 2.0);
  LC_ASSERT("Testing Lua table", sub.getNumber("c") == 3.0);

  LC_ASSERT("Testing Lua table", back.getNumber("r") == 0.3);
  LC_ASSERT("Testing Lua table", back.getNumber("b") == 0.1);
  LC_ASSERT("Testing Lua table", back.getNumber("g") == 0.0);
  LC_ASSERT("Testing Lua table", back.getString("name") == "blue_green");

  LC_RAISES("Testing Lua table", back.getNumber("name"), Lucee::Except);

  LC_ASSERT("Testing Lua table", back.hasNumVec("cells"));
  std::vector<double> cells = back.getNumVec("cells");
  LC_ASSERT("Testing list of numbers", cells.size() == 2);
  LC_ASSERT("Testing list of numbers", cells[0] == 100);
  LC_ASSERT("Testing list of numbers", cells[1] == 50);

  LC_ASSERT("Testing Lua table", back.hasStrVec("address"));
  std::vector<std::string> address = back.getStrVec("address");
  LC_ASSERT("Testing list of strings", address.size() == 2);
  LC_ASSERT("Testing list of strings", address[0] == "hello");
  LC_ASSERT("Testing list of strings", address[1] == "world");
}

void
test_2(Lucee::LuaState& L)
{
// string with table
  std::string tblStr = "coeffs = {1.0, 0.0, 0.5}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "coeffs");

// construct LuaTable object
  Lucee::LuaTable coeffs(L, "coeffs");
// test it
  std::vector<double> vals = coeffs.getAllNumbers();
  LC_ASSERT("Testing if vals has proper size", vals.size() == 3);
  LC_ASSERT("Testing if vals is correct", vals[0] == 1.0);
  LC_ASSERT("Testing if vals is correct", vals[1] == 0.0);
  LC_ASSERT("Testing if vals is correct", vals[2] == 0.5);
}

void
test_3(Lucee::LuaState& L)
{
// string with table
  std::string tblStr = 
    "background = {"
    "__kind = \"color\","
    "__type = \"Lucee\","
    "r = 0.3, b = 0.1, g = 0.0,"
    "name = \"blue_green\","
    "cells = {100, 50},"
    "address = {\"hello\", \"world\"},"
    "subcolors = {a=1.0, b=2.0, c=3.0},"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// construct LuaTable object
  Lucee::LuaTable back(L, "background");

  LC_ASSERT("Testing Lua table", back.hasNumVec("cells-nothere") == false);
  LC_ASSERT("Testing Lua table", back.hasStrVec("cells-nothere") == false);
}

void
test_4(Lucee::LuaState& L)
{
// string with table
  std::string tblStr = 
    "background = {"
    "  grida = {__type = \"Grid\"}, "
    "  gridb = {__type = \"Grid\"}, "
    "  gridc = {__type = \"Grid\"}, "
    "  solvera = {__type = \"Updater\"}, "
    "  solverb = {__type = \"Updater\"}, "
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// construct LuaTable object
  Lucee::LuaTable back(L, "background");

  std::vector<std::string> grids = back.getNamesOfType("Grid");
  LC_ASSERT("Testing if objects of Grid type are correct",
    grids.size() == 3);
  LC_ASSERT("Testing names of Grid", grids[0] == "gridc");
  LC_ASSERT("Testing names of Grid", grids[1] == "gridb");
  LC_ASSERT("Testing names of Grid", grids[2] == "grida");

  Lucee::LuaTable grida = back.getTable("grida");
  LC_ASSERT("Testing type of table", grida.getString("__type") == "Grid");

  std::vector<std::string> updaters = back.getNamesOfType("Updater");
  LC_ASSERT("Testing if objects of Updater type are correct",
    updaters.size() == 2);
  LC_ASSERT("Testing names of Updater", updaters[0] == "solvera");
  LC_ASSERT("Testing names of Updater", updaters[1] == "solverb");

  LC_ASSERT("Testing if incorrect type can be queried",
    back.getNamesOfType("not-there").size() == 0);
}

void
test_5(Lucee::LuaState& L)
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
  Lucee::LuaFunction fun(L, "init");

  double t, xyz[3];

  t = 1.0;
  xyz[0] = 1.0; xyz[1] = 2.0; xyz[2] = 3.0;
  LC_ASSERT("Testing function evaluation", fun.eval(t, xyz) == t*(1+2+3));

  t = 2.0;
  xyz[0] = 1.0; xyz[1] = 2.0; xyz[2] = 3.0;
  LC_ASSERT("Testing function evaluation", fun.eval(t, xyz) == t*(1+2+3));

  t = 3.0;
  xyz[0] = 1.0; xyz[1] = 2.0; xyz[2] = 3.0;
  LC_ASSERT("Testing function evaluation", fun.eval(t, xyz) == t*(1+2+3));

  t = 10.0;
  double dx=0.1, dy=0.2, dz=0.3;
  for (unsigned i=0; i<100; ++i)
    for (unsigned j=0; j<100; ++j)
      for (unsigned k=0; k<100; ++k)
      {
        double x = 1.0 + dx*i;
        double y = 2.0 + dy*j;
        double z = 3.0 + dz*k;
        xyz[0] = x; xyz[1] = y; xyz[2] = z;
        LC_ASSERT("Testing function evaluation", 
          fun.eval(t, xyz) == t*(x+y+z));
      }

  t = 10.0;
  for (unsigned i=0; i<100; ++i)
    for (unsigned j=0; j<100; ++j)
      {
        double x = 1.0 + dx*i;
        double y = 2.0 + dy*j;
        double z = 0.0;
        xyz[0] = x; xyz[1] = y; xyz[2] = z;
        LC_ASSERT("Testing function evaluation", 
          fun.eval(t, xyz) == t*(x+y+z));
      }
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluatable");
  Lucee::LuaState L;
  test_1(L);
  test_2(L);
  test_3(L);
  test_4(L);
  test_5(L);

  LC_END_TESTS;
}
