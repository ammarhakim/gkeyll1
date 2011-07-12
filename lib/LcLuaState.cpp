/**
 * @file	LcLuaState.cpp
 *
 * @brief	Class to represent Lua state
 */

// lucee includes
#include <LcLuaMathLib.h>
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
    {
// load all standard libraries
      luaL_openlibs(state);
// load lucee specific libraries
      Lucee::registerLuaMathLib(state);
    }
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
