/**
 * @file	LcFieldFunctionBoundaryCondition.cpp
 *
 * @brief	Class for applying boundary conditions specified via a Lua function.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFieldFunctionBoundaryCondition.h>
#include <LcGlobals.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set costructor name
  const char *FieldFunctionBoundaryCondition::id = "FieldFunction";

  void
  FieldFunctionBoundaryCondition::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);

// get list of components to apply this boundary condition
    std::vector<double> cd = tbl.getNumVec("inpComponents");
    for (unsigned i=0; i<cd.size(); ++i)
      inpComponents.push_back( (int) cd[i] );

// get reference to function
    fnRef = tbl.getFunctionRef("bc");
  }

  void
  FieldFunctionBoundaryCondition::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin1,
    const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(*L, loc[i]);
    lua_pushnumber(*L, tm);
// push
    for (unsigned i=0; i<inpComponents.size(); ++i)
      lua_pushnumber(*L, qin[inpComponents[i]]);

    unsigned numInp = 4+inpComponents.size();
    unsigned numOut = this->numComponents();
// call function
    if (lua_pcall(*L, numInp, numOut, 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("FieldFunctionBoundaryCondition::applyBc: ");
      lce << "Problem evaluating function supplied as 'bc' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-this->numComponents(); i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("FieldFunctionBoundaryCondition::applyBc: Return value not a number");
      qbc[this->component(this->numComponents()+i)] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
  }
}
