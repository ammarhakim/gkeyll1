/**
 * @file	LcLuaSaveRstrStack.cpp
 *
 * @brief	Class to save/restore LUA stack.
 *
 * @version	$Id: LcLuaSaveRstrStack.cpp 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLuaSaveRstrStack.h>

// lua includes
#include <lua.hpp>

namespace Lucee
{
  LuaSaveRstrStack::LuaSaveRstrStack(Lucee::LuaState& lcState)
    : state(lcState)
  {
    if (lcState.isValid())
      stackTop = lua_gettop(state);
  }

  LuaSaveRstrStack::~LuaSaveRstrStack()
  {
    lua_settop(state, stackTop);
  }
}
