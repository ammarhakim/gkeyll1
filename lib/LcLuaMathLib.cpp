/**
 * @file	LcLuaMathLib.cpp
 *
 * @brief	Lua callable math functions and physical constants.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcLuaMathLib.h>
#include <LcMathLib.h>
#include <LcMathPhysConstants.h>

static
int
lua_lengendre(lua_State *L)
{
  int l = (int) luaL_checknumber(L, 1);
  int m = (int) luaL_checknumber(L, 2);
  double x = luaL_checknumber(L, 3);
  lua_pushnumber(L, Lucee::legendre(l, m, x));
  return 1;
}

static const luaL_Reg lcLuaMathLib[] = {
    {"legendre", lua_lengendre},
    {NULL, NULL}
};

namespace Lucee
{
  void registerLuaMathLib(lua_State *L)
  {
// register functions
    luaL_register(L, "Lucee", lcLuaMathLib);

// register math/phys constants
    lua_pushnumber(L, Lucee::PI);
    lua_setfield(L, -2, "Pi");

    lua_pushnumber(L, Lucee::E);
    lua_setfield(L, -2, "E");

    lua_pushnumber(L, Lucee::SPEED_OF_LIGHT);
    lua_setfield(L, -2, "SpeedOfLight");

    lua_pushnumber(L, Lucee::ELECTRON_MASS);
    lua_setfield(L, -2, "ElectronMass");

    lua_pushnumber(L, Lucee::PLANCKS_CONSTANT_H);
    lua_setfield(L, -2, "PlancksConstant");
  }
}
