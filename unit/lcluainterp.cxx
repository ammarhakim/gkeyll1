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
#include <vector>

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

void
putIntoTable(lua_State *L, const std::string& key, const std::string& val)
{
  lua_pushstring(L, key.c_str());
  lua_pushstring(L, val.c_str());
  lua_settable(L, -3);
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
  if (! lua_istable(L, 1) )
    throw Lucee::Except("Function expects table as single parameter");
// add meta-data keys into table
  putIntoTable(L, "__kind", "RteHomogenousSlab");
  putIntoTable(L, "__type", "Solvers");

  return 1;
}

int
luaopen_mylib(lua_State *L)
{
  std::vector<luaL_Reg> libList;

  luaL_Reg double_lr = {"twice", l_double};
  libList.push_back(double_lr);

  luaL_Reg tableCtor_lr = {"tableCtor", l_tblctor};
  libList.push_back(tableCtor_lr);

  luaL_Reg null_lr = {NULL, NULL};
  libList.push_back(null_lr);

  luaL_register(L, "lucee", &libList[0]);
  return 1;
}

class LuaTable
{
  public:
    LuaTable(Lucee::LuaState& L)
      : ls(L) 
    {
      tblRef = luaL_ref(L, LUA_REGISTRYINDEX);
    }

    ~LuaTable()
    {
      luaL_unref(ls, LUA_REGISTRYINDEX, tblRef);
    }

    double getDouble(const std::string& key)
    {
      lua_rawgeti(ls, LUA_REGISTRYINDEX, tblRef);
      return getFromTable(ls, key);
      lua_pop(ls, 1);
    }

    LuaTable getTable(const std::string& nm)
    {
      lua_rawgeti(ls, LUA_REGISTRYINDEX, tblRef);
      lua_pushstring(ls, nm.c_str());
      lua_gettable(ls, -2); // -1 is key, -2 is table
      if (! lua_istable(ls, -1) )
      {
        Lucee::Except lce("Key ");
        lce << nm << " is not a table" << std::endl;
      }
      return LuaTable(ls);
    }

  private:
    Lucee::LuaState& ls;
    int tblRef;
};

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
  
  LuaTable ltbl(L);
  
  LC_ASSERT("Testing if table got properly", ltbl.getDouble("r") == 0.3);
  LC_ASSERT("Testing if table got properly", ltbl.getDouble("b") == 0.1);
  LC_ASSERT("Testing if table got properly", ltbl.getDouble("g") == 0.0);

  LuaTable ltbl2 = ltbl.getTable("sub_table");
  LC_ASSERT("Testing if sub-table got properly", ltbl2.getDouble("s") == 2.0);
  LC_ASSERT("Testing if sub-table got properly", ltbl2.getDouble("v") == 1.0);

  LC_ASSERT("Testing if table got properly", ltbl.getDouble("h") == 3.0);

  LC_ASSERT("Testing if function called worked", luaFunc(L, "addTwo", 2.0, 1.5) == 3.5);
  LC_ASSERT("Testing if function called worked", luaFunc(L, "minusTwo", 2.0, 1.5) == 0.5);

  LuaTable ltbl3 = ltbl2.getTable("subsub_table");
  LC_ASSERT("Testing if function called worked", ltbl3.getDouble("d") == 0.0);
  LC_ASSERT("Testing if function called worked", ltbl3.getDouble("e") == 1.0);
}

int
main(void)
{
  LC_BEGIN_TESTS("lcluainterp");
  Lucee::LuaState luaState;
  load(luaState, "lcluainterp-1.lua");

  LC_END_TESTS;
}
