/**
 * @file	LcLuaMathLib.h
 *
 * @brief	Lua callable math functions and physical constants.
 */

#ifndef LC_LUA_MATH_LIB_H
#define LC_LUA_MATH_LIB_H

// lua includes
#include <lua.hpp>

namespace Lucee
{
  void registerLuaMathLib(lua_State *L);
}

#endif // LC_LUA_MATH_LIB_H
