/**
 * @file	LcLuaState.h
 *
 * @brief	Class to represent Lua state
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
 * Opens a new instance of a Lua interpreter and closes it on
 * destruction.
 */
  class LuaState
  {
    public:
/** Create a new Lua state */
      LuaState();

/**
 * Create a new Lua state from given state.
 *
 * @param L lua state to use
 */
      LuaState(lua_State *L);

/** Destroy state */
      virtual ~LuaState();

/**
 * Get a pointer to native Lua state object.
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
/** Pointer to Lua state */
      lua_State *state;
/** Flag to indicate if we own this state */
      bool isOwner;
/** Flag to indicate if state is valid */
      bool isValidState;
  };
}

#endif // LC_LUA_STATE_H
