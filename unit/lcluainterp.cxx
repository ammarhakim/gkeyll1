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
#include <LcTest.h>

// lua include
#include <lua.hpp>

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  
  lua_State *L = lua_open();
  luaL_openlibs(L);

  char buff[256];
  while (fgets(buff, sizeof(buff), stdin) != NULL) 
  {
    int error = luaL_loadbuffer(L, buff, strlen(buff), "line") ||
      lua_pcall(L, 0, 0, 0);
    if (error) 
    {
      fprintf(stderr, "%s", lua_tostring(L, -1));
      lua_pop(L, 1);  /* pop error message from the stack */
    }
  }
  lua_close(L);

  LC_END_TESTS;
}
