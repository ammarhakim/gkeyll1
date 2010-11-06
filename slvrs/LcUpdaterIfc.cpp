/**
 * @file	LcUpdaterIfc.cpp
 *
 * @brief	Base class for updaters in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcSolverAssembly.h>
#include <LcUpdaterIfc.h>

// std includes
#include <limits>

namespace Lucee
{
// set module name
  const char *UpdaterIfc::id = "Updater";

  UpdaterIfc::UpdaterIfc()
    : Lucee::BasicObj(UpdaterIfc::id), grid(0)
  {
  }

  UpdaterIfc::~UpdaterIfc()
  {
  }

  void
  UpdaterIfc::readInput(Lucee::LuaTable& tbl)
  {
// read in grid if specified (otherwise assume it will be set by setGrid)
    if (tbl.hasObject<Lucee::GridIfc>("onGrid"))
      grid = &tbl.getObjectAsBase<Lucee::GridIfc>("onGrid");
  }

  void
  UpdaterIfc::initialize()
  {
  }

  void
  UpdaterIfc::setCurrTime(double tm) 
  {
    currTime = tm;
  }

  double
  UpdaterIfc::getCurrTime() const 
  {
    return currTime;
  }

  void
  UpdaterIfc::setGrid(const Lucee::GridIfc& grd)
  {
    grid = &grd;
  }

  void
  UpdaterIfc::setInpVars(const std::vector<const Lucee::DataStructIfc*>& dsl)
  {
    inpVars.resize(dsl.size());
     for (unsigned i=0; i<dsl.size(); ++i)
       inpVars[i] = dsl[i];
  }

  void
  UpdaterIfc::setOutVars(const std::vector<Lucee::DataStructIfc*>& dsl)
  {
    outVars.resize(dsl.size());
    for (unsigned i=0; i<dsl.size(); ++i)
      outVars[i] = dsl[i];
  }

  void
  UpdaterIfc::appendInpVarType(const std::type_info& type)
  {
    inpVarTypes.varTypes.push_back(&type);
  }

  void
  UpdaterIfc::appendOutVarType(const std::type_info& type)
  {
    outVarTypes.varTypes.push_back(&type);
  }

  void
  UpdaterIfc::setLastInpVarType(const std::type_info& type)
  {
    inpVarTypes.lastVarType = &type;
  }

  void
  UpdaterIfc::setLastOutVarType(const std::type_info& type)
  {
    outVarTypes.lastVarType = &type;
  }

  void
  UpdaterIfc::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("setCurrTime", luaSetCurrTime);
    lfm.appendFunc("advance", luaAdvance);
    lfm.appendFunc("initialize", luaInitialize);
    lfm.appendFunc("setIn", luaSetInpVars);
    lfm.appendFunc("setOut", luaSetOutVars);
  }

  int
  UpdaterIfc::luaInitialize(lua_State *L)
  {
    UpdaterIfc *updater
      = Lucee::PointerHolder<UpdaterIfc>::getObj(L);
    updater->initialize();
    return 0;
  }

  int
  UpdaterIfc::luaSetCurrTime(lua_State *L)
  {
    UpdaterIfc *updater
      = Lucee::PointerHolder<UpdaterIfc>::getObj(L);
    double t = lua_tonumber(L, 2); // current time to set
    updater->setCurrTime(t);

    return 0;
  }

  int
  UpdaterIfc::luaAdvance(lua_State *L)
  {
    UpdaterIfc *updater
      = Lucee::PointerHolder<UpdaterIfc>::getObj(L);
    double t = lua_tonumber(L, 2); // time to advance to
    Lucee::UpdaterStatus s = updater->update(t);
// push status and suggested time-step on stack
    if (s.status)
      lua_pushboolean(L, 1);
    else
      lua_pushboolean(L, 0);
    lua_pushnumber(L, s.dt);

    return 2;
  }

  int
  UpdaterIfc::luaSetInpVars(lua_State *L)
  {
    UpdaterIfc *updater
      = Lucee::PointerHolder<UpdaterIfc>::getObj(L);
    if (lua_type(L, 2) != LUA_TTABLE)
    {
      Lucee::Except lce("UpdaterIfc::luaSetInpVars: Must provide a table of datastructures to 'setInp' method");
      throw lce;
    }
    lua_pushvalue(L, 2); // push table on stack
    Lucee::LuaState ls(L);
    LuaTable tbl = Lucee::LuaTable(ls, "setInp");
    lua_pop(L, 2);
// get list of all inputs
    std::vector<const DataStructIfc*> inp = tbl.getAllObjects<const DataStructIfc>();
    updater->setInpVars(inp);
    return 0;
  }

  int
  UpdaterIfc::luaSetOutVars(lua_State *L)
  {
    UpdaterIfc *updater
      = Lucee::PointerHolder<UpdaterIfc>::getObj(L);
    if (lua_type(L, 2) != LUA_TTABLE)
    {
      Lucee::Except lce("UpdaterIfc::luaSetOutVars: Must provide a table of datastructures to 'setOut' method");
      throw lce;
    }
    lua_pushvalue(L, 2); // push table on stack
    Lucee::LuaState ls(L);
    LuaTable tbl = Lucee::LuaTable(ls, "setOut");
    lua_pop(L, 2);
// get list of all inputs
    std::vector<DataStructIfc*> out = tbl.getAllObjects<DataStructIfc>();
    updater->setOutVars(out);
    return 0;
  }
}
