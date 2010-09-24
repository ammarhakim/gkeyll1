/**
 * @file	LcLuaFuncMap.cpp
 *
 * @brief	Class to store map of callable Lua functions.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaFuncMap.h>

namespace Lucee
{
  void
  LuaFuncMap::appendFunc(const std::string& nm, int (*func)(lua_State *L))
  {
    std::map<std::string, int (*)(lua_State *L)>::iterator itr =
      funcs.find(nm);
    if (itr != funcs.end())
// replace by supplied function if already exists
      itr->second = func;
    else
    { // add a new entry to map
      funcs.insert(std::pair<std::string, int (*)(lua_State *L)>(
          nm, func));
    }
  }
}
