/**
 * @file	LcLuaFunction.cpp
 *
 * @brief	Class to represent Lua function.
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
#include <LcExcept.h>
#include <LcLuaFunction.h>

// std include
#include <iostream>

/*
#define SHOW_LUA_STACK_SIZE(nm, L)               \
    do { \
      std::cout << "Stack size in " << nm << " is " << lua_gettop(L) << std::endl; \
    } while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do { \
    std::cout << "End stack size is " << lua_gettop(L) << std::endl; \
    } while (0);
*/

#define SHOW_LUA_STACK_SIZE(nm, L)               \
    do {} while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do {} while (0);

namespace Lucee
{
  LuaFunction::LuaFunction(Lucee::LuaState& lin, const std::string& nm)
    : L(lin), name(nm)
  {
// test top of stack
    if (! lua_isfunction(L, -1) )
    {
      Lucee::Except lce("LuaFunction::LuaFunction: ");
      lce << nm << " is not a Lua function";
      throw lce;
    }
// get reference to function object
    ref = luaL_ref(L, LUA_REGISTRYINDEX);
// create typeMap
  }

  LuaFunction::~LuaFunction()
  {
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }

  double
  LuaFunction::eval(double t, double xyz[3])
  {
    SHOW_LUA_STACK_SIZE("eval", L)
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push variables on stack
    lua_pushnumber(L, t);
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(L, xyz[i]);
    if (lua_pcall(L, 4, 1, 0) != 0)
    {
      Lucee::Except lce("LuaFunction::eval: Problem evaluating function ");
      lce << name;
      throw lce;
    }
    if (!lua_isnumber(L, -1))
      throw Lucee::Except("LuaFunction::eval: Return value not a number");
    double res = lua_tonumber(L, -1);
    lua_pop(L, 1);
    SHOW_LUA_STACK_SIZE2(L)
    return res;
  }
}
