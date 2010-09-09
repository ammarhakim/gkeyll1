/**
 * @file	LcLuaState.cpp
 *
 * @brief	Class to represent Lua state
 *
 * @version	$Id: LcLuaState.cpp 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcLuaState.h>

namespace Lucee
{
  LuaState::LuaState()
    : state(lua_open()), isOwner(true)
  {
    isValidState = false;
    if (state) 
      isValidState = true;
      
    if (isValid())
// load all standard libraries
      luaL_openlibs(state);
  }
  
  LuaState::LuaState(lua_State *L)
    : state(L), isOwner(false)
  {
    isValidState = false;
    if (state) 
      isValidState = true;
  }

  LuaState::~LuaState()
  {
    if (isOwner)
    {
      if (isValid())
        lua_close(state);
    }
  }
}
