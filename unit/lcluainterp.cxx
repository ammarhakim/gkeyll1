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

double
getFromTable(Lucee::LuaState& L, const std::string& key)
{
  lua_pushstring(L, key.c_str());
  lua_gettable(L, -2); // -1 is key, -2 is table
  double res = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return res;
}

double
luaFunc(Lucee::LuaState& L, const std::string& fn, double x, double y)
{
  lua_getglobal(L, fn.c_str());
  lua_pushnumber(L, x);
  lua_pushnumber(L, y);

  if (lua_pcall(L, 2, 1, 0) != 0)
    throw Lucee::Except("Unable to call function");

  double z = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return z;
}

class MyClass
{
  public:
    static int build(lua_State *L)
    {
      std::cout << "MyClass::build" << std::endl;
      lua_pushinteger(L, 0);
      return 1;
    }

  private:

};

int
l_double(lua_State *L)
{
  double d = lua_tonumber(L, 1);
  lua_pushnumber(L, 2*d);
  return 1;
}

int
l_tblctor(lua_State *L)
{
  std::cout << "l_tblctor" << std::endl;

  lua_pushinteger(L, 0);
  return 1;
}

static const struct luaL_Reg mylib[] = 
{
  {"twice", l_double},
  {"myClass", MyClass::build},
  {"tblctor", l_tblctor},
  {NULL, NULL}
};

int
luaopen_mylib(lua_State *L)
{
  luaL_register(L, "lucee", mylib);
  return 1;
}

void
load(Lucee::LuaState& L, const std::string& fname)
{
  luaopen_mylib(L);

  if (luaL_loadfile(L, fname.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    Lucee::Except lce("Cannot read config file ");
    lce << fname << std::endl;
    throw lce;
  }
  lua_getglobal(L, "width");
  LC_ASSERT("Testing if reading from LUA worked", lua_tonumber(L, -1) == 200);
  lua_getglobal(L, "height");
  LC_ASSERT("Testing if reading from LUA worked", lua_tonumber(L, -1) == 300);
  lua_getglobal(L, "background");
  if (! lua_istable(L, -1))
    throw Lucee::Except("Table 'background' not found");
  
  LC_ASSERT("Testing if table got properly", getFromTable(L, "r") == 0.3);
  LC_ASSERT("Testing if table got properly", getFromTable(L, "b") == 0.1);
  LC_ASSERT("Testing if table got properly", getFromTable(L, "g") == 0.0);  

//   LC_ASSERT("Testing if function called worked", luaFunc(L, "addTwo", 2.0, 1.5) == 3.5);
//   LC_ASSERT("Testing if function called worked", luaFunc(L, "minusTwo", 2.0, 1.5) == 0.5);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  Lucee::LuaState luaState;
  load(luaState, "lcluainterp-1.lua");

  LC_END_TESTS;
}
