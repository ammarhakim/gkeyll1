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
#include <vector>
#include <string>

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
      std::vector<luaL_Reg> regObjFuncs;
  };
}

#endif // LC_LUA_MODULE_H
