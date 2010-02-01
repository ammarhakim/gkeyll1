/**
 * @file	LcLuaState.h
 *
 * @brief	Class to represent LUA state
 *
 * @version	$Id: LcLuaState.h 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_STATE_H
#define LC_LUA_STATE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lua includes
#include <lua.hpp>

namespace Lucee
{
/** 
 * Opens a new instance of a LUA interpreter and closes it on
 * destruction.
 */
  class LuaState
  {
    public:
/** Create a new LUA state */
      LuaState();

/** Destroy state */
      virtual ~LuaState();

/**
 * Get a pointer to native LUA state object.
 *
 * @return lua_state object.
 */
      lua_State* getState() 
      {
        return state;
      }

/**
 * Cast operator so this can be passed to functions expecting
 * lua_State*.
 *
 * @return lua_State pointer.
 */
      operator lua_State* ()
      {
        return state;
      }

/**
 * Check if state is valid.
 *
 * @return true if state is valid, false otherwise. 
 */
      bool isValid() const
      {
        return isValidState;
      }

    private:
/** Pointer to LUA state */
      lua_State *state;
/** Flag to indicate if state is valid */
      bool isValidState;
  };
}


#endif // LC_LUA_STATE_H
