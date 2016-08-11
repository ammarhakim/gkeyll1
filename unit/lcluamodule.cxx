/**
 * @file	lcluamodule.cxx
 *
 * @brief	Unit tests for using LUA into Lucee
 */

// lucee includes
#include <LcBasicObj.h>
#include <LcExcept.h>
#include <LcLuaModuleRegistry.h>
#include <LcLuaState.h>
#include <LcLuaTable.h>
#include <LcObjRegistry.h>
#include <LcTest.h>

// std includes
#include <cstdlib>

class Solver : public Lucee::BasicObj
{
  public:
    static const char *id;
    virtual ~Solver() {}

    virtual std::string what() = 0;

    virtual void readInput(Lucee::LuaTable& tbl)
    {
    }
};
const char *Solver::id = "Solver";

class RteSolver : public Solver
{
  public:
    static const char *id;

    virtual ~RteSolver() {}

    virtual std::string what()
    {
      return "RteSolver";
    }
};
const char *RteSolver::id = "RteSolver";

void registerEverything()
{
// register object
  new Lucee::ObjRegistry<Solver, RteSolver>;
}

void
test_1(Lucee::LuaState& L)
{
// string with table
  std::string tblStr = 
    "simulation = Solver.RteSolver {"
    "cells = {100, 50},"
    "}";
// evaluate string as Lua code
  if (luaL_loadstring(L, tblStr.c_str()) || lua_pcall(L, 0, 0, 0))
  {
    std::string err(lua_tostring(L, -1));
    lua_pop(L, 1);
    std::cerr << err << std::endl;
    exit(1);
  }

// put table on top of stack
  lua_getglobal(L, "simulation");
// create table
//   Solver *slvr = Lucee::ObjCreator<Solver>::getNew("RteSolver");
//   LC_ASSERT("Testing created object", slvr->what() == "RteSolver");

}

int
main(int argc, char **argv)
{
  LC_BEGIN_TESTS("lcluamodule");
  Lucee::LuaState L;

// create top-level logger
  Lucee::Logger& logger = Lucee::Logger::create("lucee");
  logger.setLevel("debug"); // base logger should log everything
// create console logger
  Lucee::Logger& conLogger = Lucee::Logger::create("lucee.console");

  registerEverything();
  //Lucee::ObjCreator<Solver>::registerModule(L);
  Lucee::LuaModuleRegistry<Solver>::registerModule(L);

  test_1(L);

  LC_END_TESTS;
}
