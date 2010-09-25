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
#include <vector>

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

/**
 * Set the deletion function to delete the object.
 *
 * @param func Pointer to Lua callable function.
 */
      void setDelFunc(int (*func)(lua_State *L));

/**
 * Get the deletion function to delete the object.
 *
 * @return Pointer to Lua callable function.
 */
      int (*getDelFunc())(lua_State *)
      {
        return delFunc;
      }

/**
 * Get list of methods added to this class in a std::vector that can
 * then be passed to Lua registration system.
 *
 * @param funcLst On output filled with Lua functions added to this class.
 */
      void fillWithFuncList(std::vector<luaL_Reg>& funcLst);

    private:
/** Map of function name to function pointer */
      std::map<std::string, int (*)(lua_State *L)> funcs;
/** Function pointer to deletion function */
      int (*delFunc)(lua_State *L);
  };
}

#endif // LC_LUA_FUNC_MAP_H
