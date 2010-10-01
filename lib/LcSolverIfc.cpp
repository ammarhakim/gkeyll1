/**
 * @file	LcSolverIfc.cpp
 *
 * @brief	Interface class for Lucee solvers.
 *
 * @version	$Id: LcSolverIfc.cpp 312 2010-03-02 18:37:24Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSolverIfc.h>

// std includes
#include <string>

namespace Lucee
{
// set madule name
  const char *SolverIfc::id = "Solver";

  SolverIfc::SolverIfc(const std::string& nm)
    : Lucee::BasicObj(nm)
  {
  }

  SolverIfc::~SolverIfc()
  {
  }

  void
  SolverIfc::setCurrTime(double tm) 
  {
    currTime = tm;
  }

  void
  SolverIfc::setNumOutFrames(unsigned n)
  {
    numOutFrames = n;
  }

  void
  SolverIfc::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  SolverIfc::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("advance", luaAdvance);
    lfm.appendFunc("initialize", luaInitialize);
    lfm.appendFunc("write", luaWrite);
  }

  int
  SolverIfc::luaAdvance(lua_State *L)
  {
    SolverIfc *solver
      = Lucee::PointerHolder<SolverIfc>::checkUserType(L);
    double t = lua_tonumber(L, 2); // time to advance to
    int res = solver->advance(t);
    lua_pushnumber(L, res);
    return 1;
  }

  int
  SolverIfc::luaInitialize(lua_State *L)
  {
    SolverIfc *solver
      = Lucee::PointerHolder<SolverIfc>::checkUserType(L);
    solver->initialize();
    return 0;
  }

  int
  SolverIfc::luaWrite(lua_State *L)
  {
    SolverIfc *solver
      = Lucee::PointerHolder<SolverIfc>::checkUserType(L);
    std::string nm = lua_tostring(L, 2); // base name
    int d = lua_tonumber(L, 3); // dump-number
    solver->writeToFile(nm, d);
    return 0;
  }
}
