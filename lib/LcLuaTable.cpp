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

#include <iostream>

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
// create typeMap
    createTypeMap();
  }

  LuaTable::~LuaTable()
  {
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }

  std::string
  LuaTable::getType()
  {
    return getString("__type");
  }

  std::string
  LuaTable::getKind()
  {
    return getString("__kind");
  }

  std::vector<double>
  LuaTable::getAllNumbers()
  {
    std::vector<double> res;
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) == LUA_TNUMBER)
        res.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
    return res;
  }

  std::string
  LuaTable::getString(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TSTRING)
    {
      Lucee::Except lce("LuaTable::getString: ");
      lce << key << " is not a string";
      throw lce;
    }
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
    if (lua_type(L, -1) != LUA_TNUMBER)
    {
      Lucee::Except lce("LuaTable::getNumber: ");
      lce << key << " is not a number";
      throw lce;
    }
    double res = lua_tonumber(L, -1);

    lua_pop(L, 1);
    return res;
  }

  std::vector<std::string>
  LuaTable::getStrVec(const std::string& key)
  {
    std::vector<std::string> res;
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);
    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TSTRING)
      {
        Lucee::Except lce("LuaTable::getStrVec: ");
        lce << key << " is not a table of strings";
        throw lce;
      }
      res.push_back(lua_tostring(L, -1));
      lua_pop(L, 1);
    }
    return res;
  }

  std::vector<double>
  LuaTable::getNumVec(const std::string& key)
  {
    std::vector<double> res;
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);
    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TNUMBER)
      {
        Lucee::Except lce("LuaTable::getNumVec: ");
        lce << key << " is not a table of numbers";
        throw lce;
      }
      res.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
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

  bool
  LuaTable::hasString(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TSTRING)
      return false;
    return true;
  }

  bool
  LuaTable::hasNumber(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNUMBER)
      return false;
    return true;
  }

  bool
  LuaTable::hasStrVec(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);

// check if this is a table in the first place
    if (lua_type(L, -1) != LUA_TTABLE)
      return false;

    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TSTRING)
        return false;
      lua_pop(L, 1);
    }
    return true;
  }

  bool
  LuaTable::hasNumVec(const std::string& key)
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);
// check if this is a table in the first place
    if (lua_type(L, -1) != LUA_TTABLE)
      return false;

    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TNUMBER)
        return false;
      lua_pop(L, 1);
    }
    return true;
  }

  bool
  LuaTable::hasTable(const std::string& nm)
  {
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
    lua_gettable(L, -2);
    if (! lua_istable(L, -1) )
      return false;
    return true;
  }

  void
  LuaTable::createTypeMap()
  {
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if ( lua_istable(L, -1) )
      {
        std::string var = lua_tostring(L, -2);
        lua_pushstring(L, "__type");
        lua_gettable(L, -2);
        if (lua_type(L, -1) == LUA_TSTRING)
        {
          addToTypeMap(var, lua_tostring(L, -1));
        }
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
  }

  std::vector<std::string>
  LuaTable::getNamesOfType(const std::string& type) const
  {
    std::map<std::string, std::vector<std::string> >::const_iterator itr
      = typeMap.find(type);
    if (itr != typeMap.end())
      return itr->second;
    Lucee::Except lce("LuaTable::getNamesOfType: Type ");
    lce << type << " does not exist in table " << name << std::endl;
    throw lce;
  }

  void
  LuaTable::addToTypeMap(const std::string& var, const std::string& type)
  {
    typeMap[type].push_back(var);
  }
}
