/**
 * @file	LcLuaHelp.h
 *
 * @brief	Class to store functions in a LUA help.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_HELP_H
#define LC_LUA_HELP_H

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
#include <iostream>
#include <map>
#include <string>

namespace Lucee
{
/**
 * Class to hold function that gets help for derived classes using the
 * "describe" function.
 */
  template <class B>
  class LuaHelp
  {
    public:
      static int luaDescribe(lua_State *L)
      {
// fetch string with name of derived class to describe
        const char *baseId = lua_tostring(L, 1);
        std::cout << "DESCRIBE: Asked to describe " << baseId << std::endl;
        return 0;
      }
  };
}

#endif // LC_LUA_HELP_H
