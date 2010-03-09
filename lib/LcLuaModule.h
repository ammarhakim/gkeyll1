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
  class LuaModule
  {
    public:
/** Create new object */
      LuaModule();

/**
 * Set the name of the module.
 *
 * @param nm Name of module.
 */
      void setName(const std::string& nm);

/**
 * Add a new function to the module.
 *
 * @param nm Name of function to add.
 * @param fptr Pointer to function.
 */
      void addFunction(const std::string& nm, int (*fptr)(lua_State *));

/**
 * Register the module into Lua. This method should be called only
 * after all funcitons are added to the module.
 *
 * @param L object representing Lua state.
 */
      void registerModule(Lucee::LuaState& L);

    private:
/** Name of module */
      std::string name;
/** Has registerModule been called */
      bool registerCalled;
/** List of functions registered in module */
      std::vector<luaL_Reg> regFuncs;
  };
}

#endif // LC_LUA_MODULE_H
