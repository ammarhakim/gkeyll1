/**
 * @file	LcLuaState.h
 *
 * @brief	Class to represent LUA state
 *
 * @version	$Id: LcLuaState.h 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
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

    private:
/** Pointer to LUA state */
      lua_State *state;
  };
}


#endif // LC_LUA_STATE_H
