/**
 * @file	lcluaobj.cxx
 *
 * @brief	Unit tests for using LUA into Lucee
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcCmdLineArgs.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>

// std includes
#include <iostream>

struct Point
{
    double x, y, z;
};

int
newPoint(lua_State *L)
{
  size_t nbytes = sizeof(Point);
  Point *p = (Point *) lua_newuserdata(L, nbytes);
  p->x = 0.0;
  p->y = 0.0;
  return 1;
}

int
setXY(lua_State *L)
{
  Point *p = (Point *) lua_touserdata(L, 1);
  double x = lua_tonumber(L, 2);
  double y = lua_tonumber(L, 3);
  p->x = x;
  p->y = y;

  return 0;
}

int
getX(lua_State *L)
{
  Point *p = (Point *) lua_touserdata(L, 1);
  lua_pushnumber(L, p->x);
  return 1;
}

int
getY(lua_State *L)
{
  Point *p = (Point *) lua_touserdata(L, 1);
  lua_pushnumber(L, p->y);
  return 1;
}

int
twice(lua_State *L)
{
  double d = lua_tonumber(L, 1);
  lua_pushnumber(L, 2*d);
  return 1;
}

struct MyFieldData
{
    double lower, upper;
    int cells;
};

void
showLuaType(int ty)
{
  switch (ty)
  {
      case LUA_TNIL:
          std::cout << "LUA_TNIL" << std::endl;
          break;
      case LUA_TBOOLEAN:
          std::cout << "LUA_TBOOLEAN" << std::endl;
          break;
      case LUA_TNUMBER:
          std::cout << "LUA_TNUMBER" << std::endl;
          break;
      case LUA_TSTRING:
          std::cout << "LUA_TSTRING" << std::endl;
          break;
      case LUA_TTABLE:
          std::cout << "LUA_TTABLE" << std::endl;
          break;
      case LUA_TTHREAD:
          std::cout << "LUA_TTHREAD" << std::endl;
          break;
      case LUA_TUSERDATA:
          std::cout << "LUA_TUSERDATA" << std::endl;
          break;
      case LUA_TFUNCTION:
          std::cout << "LUA_TFUNCTION" << std::endl;
          break;
  }
}

int
newField(lua_State *L)
{
  size_t nbytes = sizeof(MyFieldData);
  MyFieldData *fd = (MyFieldData*) lua_newuserdata(L, nbytes);
// make the table top of stack
  lua_pushvalue(L, 1);
  Lucee::LuaState myL(L);
  Lucee::LuaTable tbl(myL, "field");
  lua_pop(L, 1); // pop off table from top

  double lower = tbl.getNumber("lower");
  double upper = tbl.getNumber("upper");
  double cells = tbl.getNumber("cells");

  fd->lower = lower;
  fd->upper = upper;
  fd->cells = (int) cells;

  return 1;
}

int
getFieldLower(lua_State *L)
{
  MyFieldData *fd = (MyFieldData *) lua_touserdata(L, 1);
  lua_pushnumber(L, fd->lower);
  return 1;
}

int
getFieldUpper(lua_State *L)
{
  MyFieldData *fd = (MyFieldData *) lua_touserdata(L, 1);
  lua_pushnumber(L, fd->upper);
  return 1;
}

static const struct luaL_Reg mylib[] = {
  {"twice", twice},
  {"newPoint", newPoint},
  {"setXY", setXY},
  {"getX", getX},
  {"getY", getY},
  {"newField", newField},
  {"fieldLower", getFieldLower},
  {"fieldUpper", getFieldUpper},
  {NULL, NULL}
};

int
main(int argc, char *argv[])
{
  Lucee::CmdLineArgs cmd("lucee");
  cmd.addArg("i", "INPUT", "Input file");
// parse command line
  cmd.parse(argc, argv);

// check for input file
  if (cmd.hasArg("i") == false)
  {
    cmd.showHelp();
    exit(1);
  }
  std::string inpFile = cmd.getArg("i");

  Lucee::LuaState L;
  luaL_register(L, "lucee", mylib);

// read lua script file
  if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    std::cerr << "Error parsing input file " << inpFile << std::endl;
    exit(1);
  }

  return 0;
}
