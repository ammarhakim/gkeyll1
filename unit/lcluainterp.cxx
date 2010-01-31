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
#include <lua.hpp>

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  
  lua_State *L = lua_open();
  luaL_openlibs(L);

  lua_close(L);

  LC_END_TESTS;
}
