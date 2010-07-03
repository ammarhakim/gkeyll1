/**
 * @file	LcLuaTXYZFunction.cpp
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
#include <LcLuaTXYZFunction.h>

// std include
#include <iostream>

/*
#define SHOW_LUA_STACK_SIZE(nm, L) \
    do { \
      std::cout << "Stack size in " << nm << " is " << lua_gettop(L) << std::endl; \
    } while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do { \
    std::cout << "End stack size is " << lua_gettop(L) << std::endl; \
    } while (0);
*/

#define SHOW_LUA_STACK_SIZE(nm, L) \
    do {} while (0);

#define SHOW_LUA_STACK_SIZE2(L) \
    do {} while (0);

namespace Lucee
{
  const char *LuaTXYZFunction::id = "LuaTXYZ";

  LuaTXYZFunction::LuaTXYZFunction(Lucee::LuaState& lin, const std::string& nm, unsigned numOut)
    : Lucee::FunctionIfc(4, numOut), L(&lin), name(nm)
  {
// test top of stack
    if (! lua_isfunction(*L, -1) )
    {
      Lucee::Except lce("LuaTXYZFunction::LuaTXYZFunction: ");
      lce << nm << " is not a Lua function";
      throw lce;
    }
// get reference to function object
    ref = luaL_ref(*L, LUA_REGISTRYINDEX);
  }

  LuaTXYZFunction::LuaTXYZFunction()
    : Lucee::FunctionIfc(4, 1), L(0), name("func")
  {
  }

  LuaTXYZFunction::~LuaTXYZFunction()
  {
    luaL_unref(*L, LUA_REGISTRYINDEX, ref);
  }

  void
  LuaTXYZFunction::readInput(Lucee::LuaTable& tbl)
  {
// get hold of LuaState object
    L = &tbl.getLuaState();
    this->setInpSize(4);
    this->setOutSize(1);
    if (tbl.hasNumber("numOut"))
      this->setOutSize(tbl.getNumber("numOut"));

// get hold of reference to function
    ref = tbl.getFunctionRef("f");
  }

  std::vector<double>
  LuaTXYZFunction::eval(const std::vector<double>& inp)
  {
    SHOW_LUA_STACK_SIZE("eval", L);
    unsigned numOut = this->getOutSize();
    std::vector<double> res(numOut);
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, ref);
// push variables on stack
    for (unsigned i=0; i<inp.size(); ++i)
      lua_pushnumber(*L, inp[i]);
    if (lua_pcall(*L, inp.size(), numOut, 0) != 0)
    {
      Lucee::Except lce("LuaTXYZFunction::eval: Problem evaluating function ");
      lce << name;
      throw lce;
    }
// fetch results
    for (int i=-numOut; i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("LuaTXYZFunction::eval: Return value not a number");
      res[numOut+i] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
    SHOW_LUA_STACK_SIZE2(L)
    return res;
  }
}
