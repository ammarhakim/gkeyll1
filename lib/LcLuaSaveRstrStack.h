/**
 * @file	LcLuaSaveRstrStack.h
 *
 * @brief	Class to save/restore LUA stack.
 *
 * @version	$Id: LcLuaSaveRstrStack.h 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_SAVE_RSTR_STACK_H
#define LC_LUA_SAVE_RSTR_STACK_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaState.h>

// lua includes
#include <lua.hpp>

namespace Lucee
{
/**
 * Saves the state of the stack on creations and restores it on
 * destruction.
 */
  class LuaSaveRstrStack
  {
    public:
/**
 * Create a new object from suppiled state.
 *
 * @param lcState opened state object.
 */
      LuaSaveRstrStack(Lucee::LuaState& lcState);

/**
 * Destroy object, restoring the state back to state when object was
 * created.
 */
      virtual ~LuaSaveRstrStack();

    private:
/** Pointer to LUA state */
      lua_State *state;
/** Index to stack top */
      int stackTop;
  };
}

#endif // LC_LUA_SAVE_RSTR_STACK_H
