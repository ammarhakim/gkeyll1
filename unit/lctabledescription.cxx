/**
 * @file	lctabledescripton.cxx
 *
 * @brief	Unit tests table description object.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcTableDescription.h>
#include <LcTest.h>

void
test_1()
{
  Lucee::TableDescription eulerTbl("euler");

  double gamma;
  eulerTbl.addValue<double>("gas_gamma", 1.6)
    .setHelp("Gas adiabatic constant")
    .setMinValue(0.0)
    .setVar(&gamma);

  LC_RAISES("Testing if non-existent value can be fetched", 
    eulerTbl.getValue<double>("XXX"), Lucee::Except);

// string with table
  std::string tblStr = 
    "euler = {"
    "gas_gamma = 1.4,"
    "}";
// evaluate string as Lua code
  Lucee::LuaState L;
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "euler");

// construct LuaTable object
  Lucee::LuaTable eulerLt(L, "euler");
  eulerTbl.checkAndSet(eulerLt);
  LC_ASSERT("Checking if gamma set correctly", gamma==1.4);
}

void
test_2()
{
  Lucee::TableDescription eulerTbl("euler");

  double gamma;
  eulerTbl.addValue<double>("gas_gamma")
    .setHelp("Gas adiabatic constant")
    .setMinValue(0.0)
    .setVar(&gamma);

  std::vector<int> vals(3);
  vals[0] = 1; vals[1] = 2; vals[2] = 3;
  eulerTbl.addValue<int>("eos_type")
    .setHelp("Equation-of-state type")
    .setOneOf(vals);

  LC_RAISES("Testing if non-existent value can be fetched", 
    eulerTbl.getValue<double>("XXX"), Lucee::Except);

// string with table
  std::string tblStr = 
    "euler = {"
    "gas_gamma = -1.4,"
    "eos_type = 6,"
    "}";
// evaluate string as Lua code
  Lucee::LuaState L;
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
    throw Lucee::Except("Unable to parse Lua string");

// fetch table and put on top of stack
  lua_getglobal(L, "euler");

// construct LuaTable object
  Lucee::LuaTable eulerLt(L, "euler");
  LC_RAISES("Testing if correct exception raised", eulerTbl.checkAndSet(eulerLt), Lucee::Except);

  try
  {
    eulerTbl.checkAndSet(eulerLt);
  }
  catch (Lucee::Except& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void
test_3()
{
  Lucee::TableDescription grd("grid");
  std::vector<int> cells;
  grd.addVector<int>("cells")
    .setHelp("Number of cells in each direction")
    .setLength(3)
    .setMinValue(1)
    .setVar(&cells);

  std::vector<double> lower;
  grd.addVector<double>("lower")
    .setHelp("Coordinates of lower corner of box")
    .setLength(3)
    .setVar(&lower);

  std::vector<double> upper;
  grd.addVector<double>("upper")
    .setHelp("Coordinates of upper corner of box")
    .setLength(3)
    .setVar(&upper);
}

int
main(void)
{
  LC_BEGIN_TESTS("lctabledescription");
  test_1();
  test_2();
  test_3();
  LC_END_TESTS;
}
