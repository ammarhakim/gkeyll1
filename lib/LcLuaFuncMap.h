/**
 * @file	LcLuaFuncMap.h
 *
 * @brief	Class to store map of callable Lua functions.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_FUNC_MAP_H
#define LC_LUA_FUNC_MAP_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lua includes
#include <lua.hpp>

// std includes
#include <map>
#include <string>

namespace Lucee
{
/**
 * Class to hold map of Lua callable functions.
 */
  class LuaFuncMap
  {
    public:
/**
 * Append a method to the function map.
 *
 * @param nm Name of function.
 * @param func Pointer to Lua callable function.
 */
      void appendFunc(const std::string& nm, int (*func)(lua_State *L));

    private:
/** Map of function */
      std::map<std::string, int (*)(lua_State *L)> funcs;
  };
}

#endif // LC_LUA_FUNC_MAP_H
