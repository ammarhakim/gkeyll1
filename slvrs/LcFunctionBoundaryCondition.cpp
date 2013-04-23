/**
 * @file	LcFunctionBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFunctionBoundaryCondition.h>
#include <LcGlobals.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set costructor name
  const char *FunctionBoundaryCondition::id = "Function";

  void
  FunctionBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);
// get reference to function
    fnRef = tbl.getFunctionRef("bc");
  }

  void
  FunctionBoundaryCondition::applyBc(double tm, const double loc[3],
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(*L, loc[i]);
    lua_pushnumber(*L, tm);
// call function
    if (lua_pcall(*L, 4, components.size(), 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("FunctionBoundaryCondition::applyBc: ");
      lce << "Problem evaluating function supplied as 'bc' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-components.size(); i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("FunctionBoundaryCondition::applyBc: Return value not a number");
      qbc[components[components.size()+i]] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
  }
}
