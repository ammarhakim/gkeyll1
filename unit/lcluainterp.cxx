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
#include <LcLuaState.h>
#include <LcTest.h>

// lua include
#include <lua.hpp>

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  
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

  LC_END_TESTS;
}
