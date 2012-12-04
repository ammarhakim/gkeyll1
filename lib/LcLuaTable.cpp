/**
 * @file	LcLuaTable.cpp
 *
 * @brief	Class to represent Lua table.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaTable.h>

#include <iostream>

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

  Lucee::LuaState&
  LuaTable::getLuaState()
  {
    return L;
  }

  std::string
  LuaTable::getType() const
  {
    SHOW_LUA_STACK_SIZE("getType", L);
    return getString("__type");
  }

  std::string
  LuaTable::getKind() const
  {
    SHOW_LUA_STACK_SIZE("getKind", L);
    return getString("__kind");
  }

  std::vector<double>
  LuaTable::getAllNumbers() const
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

  std::vector<std::string>
  LuaTable::getAllStrings() const
  {
    SHOW_LUA_STACK_SIZE("getAllStrings", L);
    std::vector<std::string> res;
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    int t = lua_gettop(L);
    lua_pushnil(L);
    while (lua_next(L, t) != 0) 
    {
      if (lua_type(L, -1) == LUA_TSTRING)
        res.push_back(lua_tostring(L, -1));
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  std::string
  LuaTable::getString(const std::string& key) const
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
  LuaTable::getNumber(const std::string& key) const
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

  bool
  LuaTable::getBool(const std::string& key) const
  {
    SHOW_LUA_STACK_SIZE("getBool", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TBOOLEAN)
    {
      lua_pop(L, 2);
      Lucee::Except lce("LuaTable::getBool: ");
      lce << key << " is not a bool";
      throw lce;
    }
    int res = lua_toboolean(L, -1);

    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res == 1 ? true : false;
  }

  std::vector<std::string>
  LuaTable::getStrVec(const std::string& key) const
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
  LuaTable::getNumVec(const std::string& key) const
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

  std::vector<bool>
  LuaTable::getBoolVec(const std::string& key) const
  {
    SHOW_LUA_STACK_SIZE("getBoolVec", L);
    std::vector<bool> res;
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
      if (lua_type(L, -1) != LUA_TBOOLEAN)
      {
        lua_pop(L, 2);
        Lucee::Except lce("LuaTable::getBoolVec: ");
        lce << key << " is not a table of booleans";
        throw lce;
      }
      res.push_back(lua_toboolean(L, -1) == 1 ? true : false);
      lua_pop(L, 1);
    }
    lua_pop(L, 2);
    SHOW_LUA_STACK_SIZE2(L);
    return res;
  }

  LuaTable
  LuaTable::getTable(const std::string& nm) const
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

  int
  LuaTable::getFunctionRef(const std::string& nm) const
  {
    SHOW_LUA_STACK_SIZE("getFunctionRef", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
// get data from table
    lua_gettable(L, -2);
    if (! lua_isfunction(L, -1) )
    {
      lua_pop(L, 1);
      Lucee::Except lce("LuaTable::getFunctionRef: ");
      lce << nm << " is not a Lua function";
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_pop(L, 1);
    SHOW_LUA_STACK_SIZE2(L);
    return fnRef;
  }

  bool
  LuaTable::hasString(const std::string& key) const
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
  LuaTable::hasNumber(const std::string& key) const
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
  LuaTable::hasBool(const std::string& key) const
  {
    SHOW_LUA_STACK_SIZE("hasBool", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, key.c_str());
// get data from table
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TBOOLEAN)
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
  LuaTable::hasStrVec(const std::string& key) const
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
  LuaTable::hasNumVec(const std::string& key) const
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
  LuaTable::hasBoolVec(const std::string& key) const
  {
    SHOW_LUA_STACK_SIZE("hasBoolVec", L);
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
      if (lua_type(L, -1) != LUA_TBOOLEAN)
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
  LuaTable::hasTable(const std::string& nm) const
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

  bool
  LuaTable::hasFunction(const std::string& nm) const
  {
    SHOW_LUA_STACK_SIZE("hasFunction", L);
// push table object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
    lua_pushstring(L, nm.c_str());
// get data from table
    lua_gettable(L, -2);
    if (! lua_isfunction(L, -1) )
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
