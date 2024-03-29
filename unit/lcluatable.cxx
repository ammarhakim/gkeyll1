/**
 * @file	lcluatable.cxx
 *
 * @brief	Unit tests for using LUA into Lucee
 */

// lucee includes
#include <LcExcept.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::LuaState L;
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
    "correct = true,"
    "screw = false,"
    "boolList = {true, false, true},"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// construct LuaTable object
  Lucee::LuaTable back(L, "background");

  // LC_ASSERT("Testing kind field", back.getKind() == "color");
  // LC_ASSERT("Testing type field", back.getType() == "Lucee");

// test it
  LC_ASSERT("Testing Lua table", back.hasNumber("r"));
  LC_ASSERT("Testing Lua table", back.hasNumber("b"));
  LC_ASSERT("Testing Lua table", back.hasNumber("g"));
  LC_ASSERT("Testing Lua table", back.hasString("name"));
  LC_ASSERT("Testing Lua table", back.hasNumber("zzz") == false);
  LC_ASSERT("Testing Lua table", back.hasBool("correct"));
  LC_ASSERT("Testing Lua table", back.hasBool("screw"));
  LC_ASSERT("Testing Lua table", back.hasBool("drew") == false);

  LC_ASSERT("Testing Lua table", back.getNumber("r") == 0.3);
  LC_ASSERT("Testing Lua table", back.getNumber("b") == 0.1);
  LC_ASSERT("Testing Lua table", back.getNumber("g") == 0.0);
  LC_ASSERT("Testing Lua table", back.getString("name") == "blue_green");
  LC_ASSERT("Testing Lua table", back.getBool("correct") == true);
  LC_ASSERT("Testing Lua table", back.getBool("screw") == false);

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

  LC_ASSERT("Testing Lua table", back.hasBoolVec("boolList"));
  std::vector<bool> blist = back.getBoolVec("boolList");
  LC_ASSERT("Testing list of bools", blist.size() == 3);
  LC_ASSERT("Testing list of bools", blist[0] == true);
  LC_ASSERT("Testing list of bools", blist[1] == false);
  LC_ASSERT("Testing list of bools", blist[2] == true);

}

void
test_2()
{
  Lucee::LuaState L;
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
test_3()
{
  Lucee::LuaState L;
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
test_4()
{
  Lucee::LuaState L;
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
test_5()
{
  Lucee::LuaState L;
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
}

void
test_6()
{
  Lucee::LuaState L;
// string with table
  std::string tblStr = 
    "myNumbers = {1,2,3,4,5}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "myNumbers");

// construct LuaTable object
  Lucee::LuaTable back(L, "myNumbers");
  std::vector<double> nums = back.getAllNumbers();
  for (unsigned i=0; i<nums.size(); ++i)
    LC_ASSERT("Tesitng list of numbers", nums[i] == i+1);
}

void
test_7()
{
  Lucee::LuaState L;
// string with table
  std::string tblStr = 
    "myNumbers = {\"a\", \"b\", \"c\"}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "myNumbers");

// construct LuaTable object
  Lucee::LuaTable back(L, "myNumbers");
  std::vector<std::string> strs = back.getAllStrings();
  LC_ASSERT("Tesitng list of strings", strs[0] == "a");
  LC_ASSERT("Tesitng list of strings", strs[1] == "b");
  LC_ASSERT("Tesitng list of strings", strs[2] == "c");
}

void
test_8()
{
  Lucee::LuaState L;
// string with table
  std::string tblStr = 
    "funcTbl = {"
    "  source = function (x, y, z)"
    "    return x, y, z"
    "  end,"
    "  gas_gamma = 1.4,"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "funcTbl");

// construct LuaTable object
  Lucee::LuaTable tbl(L, "funcTbl");

  LC_ASSERT("Checking gas gamma", tbl.getNumber("gas_gamma") == 1.4);
  LC_ASSERT("Checking for function to set source", tbl.hasFunction("source"));

  double xc[3] = {1.5, 2.5, 3.5};
  double out[3] = {0.0, 0.0, 0.0};
  int fRef = tbl.getFunctionRef("source");
// push function object on stack
  lua_rawgeti(L, LUA_REGISTRYINDEX, fRef);
// push variables on stack
  for (unsigned i=0; i<3; ++i)
    lua_pushnumber(L, xc[i]);
  if (lua_pcall(L, 3, 3, 0) != 0)
  {
    Lucee::Except lce("lcluatable: ");
    lce << "Problem evaluating function supplied to 'set' method";
    throw lce;
  }
// fetch results
  for (int i=-3; i<0; ++i)
  {
    if (!lua_isnumber(L, i))
      throw Lucee::Except("lcluatable: Return value not a number");
    out[3+i] = lua_tonumber(L, i);
  }
  lua_pop(L, 1);
  LC_ASSERT("Testing function call", out[0] == 1.5);
  LC_ASSERT("Testing function call", out[1] == 2.5);
  LC_ASSERT("Testing function call", out[2] == 3.5);
}

double
callLuaFunction(Lucee::LuaState& L, int fnRef, double param)
{
// push function object on stack
  lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
  lua_pushnumber(L, param);
  if (lua_pcall(L, 1, 1, 0) != 0)
  {
    Lucee::Except lce("lcluatable: ");
    lce << "Problem evaluating function supplied to 'callLuaFunction' method";
    throw lce;
  }
// fetch results
  if (!lua_isnumber(L, -1))
    throw Lucee::Except("lcluatable: Return value not a number");
  double res = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return res;
}

void
test_9()
{
  Lucee::LuaState L;
// string with table with list of functions
  std::string tblStr = 
    "funcTbl = {"
    "  function (x) return 1*x end,"
    "  function (x) return 2*x end,"
    "  function (x) return 3*x end,"
    "  function (x) return 4*x end,"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "funcTbl");

// construct LuaTable object, and fetch list of function references
  Lucee::LuaTable tbl(L, "funcTbl");
  std::vector<int> fnRefs = tbl.getAllFunctionRefs();

  LC_ASSERT("Checking is number of functions is correct", fnRefs.size() == 4);
// test functions
  LC_ASSERT("Testing function in list", callLuaFunction(L, fnRefs[0], 2.0) == 2.0);
  LC_ASSERT("Testing function in list", callLuaFunction(L, fnRefs[1], 2.0) == 4.0);
  LC_ASSERT("Testing function in list", callLuaFunction(L, fnRefs[2], 2.0) == 6.0);
  LC_ASSERT("Testing function in list", callLuaFunction(L, fnRefs[3], 2.0) == 8.0);
}

void
test_10()
{
  Lucee::LuaState L;
// string with table
  std::string tblStr = 
    "background = { {100, 50}, {200, 100} }";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "background");

// construct LuaTable object
  Lucee::LuaTable back(L, "background");
// get first table  
  Lucee::LuaTable tbl1 = back.getTable(0);
// check its values
  std::vector<double> v = tbl1.getAllNumbers();
  LC_ASSERT("Checking length of list", v.size() == 2);
  LC_ASSERT("Checking list", v[0] == 100);
  LC_ASSERT("Checking list", v[1] == 50);
}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcluatable");

  test_1();
  test_2();
  test_3();
  test_4();
  test_5();
  test_6();
  test_7();
  test_8();
  test_9();
// test_10 is failing. Need to fix
  test_10(); 

  LC_END_TESTS;
}
