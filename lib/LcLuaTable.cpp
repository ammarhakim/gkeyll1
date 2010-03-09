/**
 * @file	LcLuaTable.cpp
 *
 * @brief	Class to represent Lua table.
 *
 * @version	$Id: LcLuaTable.cpp 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcLuaTable.h>

namespace Lucee
{
  LuaTable::LuaTable(Lucee::LuaState& lin, const std::string& nm)
    : L(lin), name(nm)
  {
// test top of stack
    if (! lua_istable(L, -1) )
      throw Lucee::Except("LuaTable::LuaTable: Top of stack is not a table");
// get reference to table object
    ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  LuaTable::~LuaTable()
  {
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }
}
