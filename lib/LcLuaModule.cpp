/**
 * @file	LcLuaModule.cpp
 *
 * @brief	Class to store functions in a LUA module.
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
#include <LcLuaModule.h>

namespace Lucee
{
  LuaModule::LuaModule()
    : registerCalled(false)
  {
  }

  void
  LuaModule::setName(const std::string& nm)
  {
    name = nm;
  }

  void
  LuaModule::addFunction(const std::string& nm, int (*fptr)(lua_State *))
  {
    luaL_Reg reg = {nm.c_str(), fptr};
    regFuncs.push_back(reg);
  }
  
  void
  LuaModule::registerModule(Lucee::LuaState& L)
  {
// do not do anything of register already called
    if (registerCalled) return;

// push back a NULL sentinel
    luaL_Reg reg = {NULL, NULL};
    regFuncs.push_back(reg);
// now register module
    luaL_register(L, name.c_str(), &regFuncs[0]);
    registerCalled = true;
  }
}
