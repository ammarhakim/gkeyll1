/**
 * @file	lcluainterp.cxx
 *
 * @brief	Unit tests for using LUA into Lucee
 *
 * @version	$Id: lcluainterp.cxx 162 2009-08-28 20:07:10Z a.hakim777 $
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcLuaState.h>
#include <LcTest.h>

// lua include
#include <lua.hpp>

// std includes
#include <string>
#include <cstring>

void
interp()
{
// create a new LUA state
  Lucee::LuaState state;

  char buff[256];
  while (fgets(buff, sizeof(buff), stdin) != NULL) 
  {
    int error = luaL_loadbuffer(state, buff, strlen(buff), "line") ||
      lua_pcall(state, 0, 0, 0);
    if (error) 
    {
      fprintf(stderr, "%s", lua_tostring(state, -1));
      lua_pop(state, 1);  /* pop error message from the stack */
    }
  }
}

void
load(Lucee::LuaState& L, const std::string& fname)
{
  if (luaL_loadfile(L, fname.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    Lucee::Except lce("Cannot read config file ");
    lce << fname << std::endl;
    throw lce;
  }
  lua_getglobal(L, "width");
  lua_getglobal(L, "height");
  LC_ASSERT("Testing if reading from LUA worked", lua_tointeger(L, -2) == 200);
  LC_ASSERT("Testing if reading from LUA worked", lua_tointeger(L, -1) == 300);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  Lucee::LuaState luaState;
  load(luaState, "lcluainterp-1.lua");

  LC_END_TESTS;
}
