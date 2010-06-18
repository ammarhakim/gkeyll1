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
  LuaFunction::LuaFunction(Lucee::LuaState& lin, const std::string& nm, unsigned numOut)
    : L(lin), name(nm), numOut(numOut)
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

  std::vector<double>
  LuaFunction::eval(const std::vector<double>& inp)
  {
    SHOW_LUA_STACK_SIZE("eval", L);

    std::vector<double> res(numOut);
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push variables on stack
    for (unsigned i=0; i<inp.size(); ++i)
      lua_pushnumber(L, inp[i]);
    if (lua_pcall(L, inp.size(), numOut, 0) != 0)
    {
      Lucee::Except lce("LuaFunction::eval: Problem evaluating function ");
      lce << name;
      throw lce;
    }
// fetch results
    for (int i=-numOut; i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("LuaFunction::eval: Return value not a number");
      res[numOut+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
    SHOW_LUA_STACK_SIZE2(L)
    return res;
  }
}
