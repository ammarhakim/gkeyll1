/**
 * @file	LcLuaModule.h
 *
 * @brief	Class to store functions in a LUA module.
 */

#ifndef LC_LUA_MODULE_H
#define LC_LUA_MODULE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaFuncMap.h>
#include <LcLuaState.h>

// lua includes
#include <lua.hpp>

// std includes
#include <map>
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Class to hold list of Lua functions that are then registered into a
 * module. These then become available from Lua script.
 */
  template <class B>
  class LuaModule
  {
    public:
/** List of registered functions to make object */
      std::vector<luaL_Reg> regCreateFuncs;
/** Map of derived class name to Lua callable functions */
      std::map<std::string, Lucee::LuaFuncMap> funcMaps;
  };
}

#endif // LC_LUA_MODULE_H
