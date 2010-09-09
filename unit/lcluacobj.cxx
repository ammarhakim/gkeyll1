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

int
newField(lua_State *L)
{
  size_t nbytes = sizeof(MyFieldData);
  MyFieldData *fd = (MyFieldData*) lua_newuserdata(L, nbytes);
    
  //Lucee::LuaState myL(L);
  //Lucee::LuaTable lt(myL, "field");
  std::cout << "Stack size in is " << lua_gettop(L) << std::endl;

  return 1;
}

int
getFieldLower(lua_State *L)
{
  MyFieldData *fd = (MyFieldData *) lua_touserdata(L, 1);
  std::cout << "Inside getFieldLower" << std::endl;
  lua_pushnumber(L, fd->lower);
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
