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


int
main(void)
{
  LC_BEGIN_TESTS("lcluatable");
  Lucee::LuaState L;
  test_1(L);
  test_2(L);

  LC_END_TESTS;
}
