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

/*
#define SHOW_LUA_STACK_SIZE(nm, L)               \
    do { \
      std::cout << "Stack size in " << nm << " is " << lua_gettop(L) << std::endl; \
    } while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do { \
    std::cout << "End stack size is " << lua_gettop(L) << std::endl; \
    } while (0);
*/

#define SHOW_LUA_STACK_SIZE(nm, L)               \
    do {} while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do {} while (0);

namespace Lucee
{
  LuaTable::LuaTable(Lucee::LuaState& lin, const std::string& nm)
    : L(lin), name(nm)
  {
    SHOW_LUA_STACK_SIZE("LuaTable", L);
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
    SHOW_LUA_STACK_SIZE2(L);
  }

  LuaTable::~LuaTable()
  {
    SHOW_LUA_STACK_SIZE("~LuaTable", L);
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
    SHOW_LUA_STACK_SIZE2(L);
  }

  std::string
  LuaTable::getType()
  {
    SHOW_LUA_STACK_SIZE("getType", L);
    return getString("__type");
  }

  std::string
  LuaTable::getKind()
  {
    SHOW_LUA_STACK_SIZE("getKind", L);
    return getString("__kind");
  }

  std::vector<double>
  LuaTable::getAllNumbers()
  {
    SHOW_LUA_STACK_SIZE("getAllNumbers", L);
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
    lua_pop(L, 1);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  std::string
  LuaTable::getString(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("getString", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TSTRING)
    {
      lua_pop(L, 2);
      Lucee::Except lce("LuaTable::getString: ");
      lce << key << " is not a string";
      throw lce;
    }
    std::string res = lua_tostring(L, -1);

    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  double
  LuaTable::getNumber(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("getNumber", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNUMBER)
    {
      lua_pop(L, 2);
      Lucee::Except lce("LuaTable::getNumber: ");
      lce << key << " is not a number";
      throw lce;
    }
    double res = lua_tonumber(L, -1);

    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  std::vector<std::string>
  LuaTable::getStrVec(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("getStrVec", L);
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
        lua_pop(L, 2);
        Lucee::Except lce("LuaTable::getStrVec: ");
        lce << key << " is not a table of strings";
        throw lce;
      }
      res.push_back(lua_tostring(L, -1));
      lua_pop(L, 1);
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  std::vector<double>
  LuaTable::getNumVec(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("getNumVec", L);
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
        lua_pop(L, 2);
        Lucee::Except lce("LuaTable::getNumVec: ");
        lce << key << " is not a table of numbers";
        throw lce;
      }
      res.push_back(lua_tonumber(L, -1));
      lua_pop(L, 1);
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  LuaTable
  LuaTable::getTable(const std::string& nm)
  {
    SHOW_LUA_STACK_SIZE("getTable", L);
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
    lua_gettable(L, -2);
    LuaTable tbl = LuaTable(L, nm);
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return tbl;
  }

  bool
  LuaTable::hasString(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("hasString", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TSTRING)
    {
      lua_pop(L, 2);
      SHOW_LUA_STACK_SIZE2(L);
      return false;
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return true;
  }

  bool
  LuaTable::hasNumber(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("hasNumber", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNUMBER)
    {
      lua_pop(L, 2);
      SHOW_LUA_STACK_SIZE2(L);
      return false;
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return true;
  }

  bool
  LuaTable::hasStrVec(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("hasStrVec", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);

// check if this is a table in the first place
    if (lua_type(L, -1) != LUA_TTABLE)
    {
      lua_pop(L, 2);
      SHOW_LUA_STACK_SIZE2(L);      
      return false;
    }

    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TSTRING)
      {
        lua_pop(L, 2);
        return false;
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return true;
  }

  bool
  LuaTable::hasNumVec(const std::string& key)
  {
    SHOW_LUA_STACK_SIZE("hasNumVec", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push name of table onto stack
    lua_pushstring(L, key.c_str());
// now put table of values on stack
    lua_gettable(L, -2);
// check if this is a table in the first place
    if (lua_type(L, -1) != LUA_TTABLE)
    {
      lua_pop(L, 2);
      SHOW_LUA_STACK_SIZE2(L);
      return false;
    }

    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) != LUA_TNUMBER)
      {
        lua_pop(L, 2);
        return false;
      }
      lua_pop(L, 1);
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return true;
  }

  bool
  LuaTable::hasTable(const std::string& nm)
  {
    SHOW_LUA_STACK_SIZE("hasTable", L);
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
    lua_gettable(L, -2);
    if (! lua_istable(L, -1) )
    {
      lua_pop(L, 2);
      SHOW_LUA_STACK_SIZE2(L);
      return false;
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return true;
  }

  void
  LuaTable::createTypeMap()
  {
//    SHOW_LUA_STACK_SIZE("createTypeMap", L);
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
//    SHOW_LUA_STACK_SIZE2(L);
  }

  std::vector<std::string>
  LuaTable::getNamesOfType(const std::string& type) const
  {
    std::map<std::string, std::vector<std::string> >::const_iterator itr
      = typeMap.find(type);
    if (itr != typeMap.end())
      return itr->second;
    return std::vector<std::string>();
  }

  void
  LuaTable::addToTypeMap(const std::string& var, const std::string& type)
  {
    typeMap[type].push_back(var);
  }
}
