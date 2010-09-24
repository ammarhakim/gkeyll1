/**
 * @file	LcLuaModule.h
 *
 * @brief	Class to store functions in a LUA module.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_MODULE_H
#define LC_LUA_MODULE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
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
/** List of registered functions */
      std::vector<luaL_Reg> regFuncs;
/** List of registered functions to make object */
      std::vector<luaL_Reg> regCreateFuncs;

/** 
 * Struct to store function name to pointer to function.
 */
      struct FuncData
      {
/** Map of function */
          std::map<std::string, int (*func)(lua_State *L)> funcs;
      };

/** Map of dervied class IDs to Lua callable functions */
      std::map<std::string, std::vector<FuncData> > callableFuncs;
  };
}

#endif // LC_LUA_MODULE_H
