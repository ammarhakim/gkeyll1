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

// test it
  LC_ASSERT("Testing Lua table", back.getNumber("r") == 0.3);
  LC_ASSERT("Testing Lua table", back.getNumber("b") == 0.1);
  LC_ASSERT("Testing Lua table", back.getNumber("g") == 0.0);
  LC_ASSERT("Testing Lua table", back.getString("name") == "blue_green");

// get subtable
  Lucee::LuaTable sub = back.getTable("subcolors");

  LC_ASSERT("Testing Lua table", sub.getNumber("a") == 1.0);
  LC_ASSERT("Testing Lua table", sub.getNumber("b") == 2.0);
  LC_ASSERT("Testing Lua table", sub.getNumber("c") == 3.0);

  LC_ASSERT("Testing Lua table", back.getNumber("r") == 0.3);
  LC_ASSERT("Testing Lua table", back.getNumber("b") == 0.1);
  LC_ASSERT("Testing Lua table", back.getNumber("g") == 0.0);
  LC_ASSERT("Testing Lua table", back.getString("name") == "blue_green");

  LC_RAISES("Testing Lua table", back.getNumber("name"), Lucee::Except);

  std::vector<double> cells = back.getNumVec("cells");
  LC_ASSERT("Testing list of numbers", cells.size() == 2);
  LC_ASSERT("Testing list of numbers", cells[0] == 100);
  LC_ASSERT("Testing list of numbers", cells[1] == 50);

  std::vector<std::string> address = back.getStrVec("address");
  LC_ASSERT("Testing list of strings", address.size() == 2);
  LC_ASSERT("Testing list of strings", address[0] == "hello");
  LC_ASSERT("Testing list of strings", address[1] == "world");
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluatable");
  Lucee::LuaState L;
  test_1(L);

  LC_END_TESTS;
}
