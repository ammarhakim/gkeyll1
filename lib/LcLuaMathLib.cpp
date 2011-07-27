/**
 * @file	LcLuaMathLib.cpp
 *
 * @brief	Lua callable math functions and physical constants.
*/

// lucee includes
#include <LcLogStream.h>
#include <LcLogger.h>
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

static
int
lua_log_debug(lua_State *L)
{
  const char *str = luaL_checkstring(L, 1);
  Lucee::Logger& l = Lucee::Logger::get("lucee.console");
  Lucee::LogStream strm = l.getDebugStream();
  strm << str << std::endl;
  return 0;
}

static
int
lua_log_info(lua_State *L)
{
  const char *str = luaL_checkstring(L, 1);
  Lucee::Logger& l = Lucee::Logger::get("lucee.console");
  Lucee::LogStream strm = l.getInfoStream();
  strm << str << std::endl;
  return 0;
}

static
int
lua_log_warning(lua_State *L)
{
  const char *str = luaL_checkstring(L, 1);
  Lucee::Logger& l = Lucee::Logger::get("lucee.console");
  Lucee::LogStream strm = l.getWarningStream();
  strm << str << std::endl;
  return 0;
}

static
int
lua_log_error(lua_State *L)
{
  const char *str = luaL_checkstring(L, 1);
  Lucee::Logger& l = Lucee::Logger::get("lucee.console");
  Lucee::LogStream strm = l.getErrorStream();
  strm << str << std::endl;
  return 0;
}

static
int
lua_log_critical(lua_State *L)
{
  const char *str = luaL_checkstring(L, 1);
  Lucee::Logger& l = Lucee::Logger::get("lucee.console");
  Lucee::LogStream strm = l.getCriticalStream();
  strm << str << std::endl;
  return 0;
}

static const luaL_Reg lcLuaMathLib[] = {
    {"legendre", lua_lengendre},
    {"logDebug", lua_log_debug},
    {"logInfo", lua_log_info},
    {"logWarning", lua_log_warning},
    {"logError", lua_log_error},
    {"logCritical", lua_log_critical},
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

    lua_pushnumber(L, Lucee::EULER);
    lua_setfield(L, -2, "Euler");

    lua_pushnumber(L, Lucee::SPEED_OF_LIGHT);
    lua_setfield(L, -2, "SpeedOfLight");

    lua_pushnumber(L, Lucee::PLANCKS_CONSTANT_H);
    lua_setfield(L, -2, "PlancksConstant");

    lua_pushnumber(L, Lucee::ELECTRON_MASS);
    lua_setfield(L, -2, "ElectronMass");

    lua_pushnumber(L, Lucee::PROTON_MASS);
    lua_setfield(L, -2, "ProtonMass");

    lua_pushnumber(L, Lucee::ELEMENTARY_CHARGE);
    lua_setfield(L, -2, "ElementaryCharge");

    lua_pushnumber(L, Lucee::BOLTZMANN_CONSTANT);
    lua_setfield(L, -2, "BoltzmannConstant");

    lua_pushnumber(L, Lucee::EPSILON0);
    lua_setfield(L, -2, "Epsilon0");

    lua_pushnumber(L, Lucee::MU0);
    lua_setfield(L, -2, "Mu0");

    lua_pop(L, 1); // pop what was pushed
  }
}
