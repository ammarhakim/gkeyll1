/**
 * @file	LcFunctionSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFunctionSource.h>
#include <LcGlobals.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for creators
  const char *FunctionSource::id = "Function";

  FunctionSource::FunctionSource()
    : Lucee::PointSourceIfc(0, 1, true)
  { 
// does not take any variables and supplies arbitrary outputs
  }

  void
  FunctionSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// get reference to function
    fnRef = tbl.getFunctionRef("source");
  }

  void
  FunctionSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(*L, loc[i]);
    lua_pushnumber(*L, tm);
// call function
    if (lua_pcall(*L, 4, src.size(), 0) != 0)
    {
      Lucee::Except lce("FunctionSource::getSource: ");
      lce << "Problem evaluating function supplied as 'source' ";
      throw lce;
    }
// fetch results
    for (int i=-src.size(); i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("FunctionSource::getSource: Return value not a number");
      src[src.size()+i] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
  }
}
