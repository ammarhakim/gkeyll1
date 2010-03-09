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
    {
      Lucee::Except lce("LuaTable::LuaTable: ");
      lce << nm << " is not a Lua table";
      throw lce;
    }
// get reference to table object
    ref = luaL_ref(L, LUA_REGISTRYINDEX);
  }

  LuaTable::~LuaTable()
  {
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }

  std::string
  LuaTable::getString(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    std::string res = lua_tostring(L, -1);

    lua_pop(L, 1);
    return res;
  }

  double
  LuaTable::getNumber(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    double res = lua_tonumber(L, -1);

    lua_pop(L, 1);
    return res;
  }

  LuaTable
  LuaTable::getTable(const std::string& nm)
  {
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
    lua_gettable(L, -2);
    return LuaTable(L, nm);
  }
}
