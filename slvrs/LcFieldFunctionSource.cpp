/**
 * @file	LcFieldFunctionSource.cpp
 *
 * @brief	Source for computing Lorentz force on a fluid
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFieldFunctionSource.h>
#include <LcGlobals.h>
#include <LcLuaState.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set id for creators
  const char *FieldFunctionSource::id = "Function";

  FieldFunctionSource::FieldFunctionSource()
    : Lucee::PointSourceIfc(1, 1, true)
  { 
// does not take any variables and supplies arbitrary outputs
  }

  void
  FieldFunctionSource::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::PointSourceIfc::readInput(tbl);
// get reference to function
    fnRef = tbl.getFunctionRef("source");
  }

  void
  FieldFunctionSource::getSource(double tm, const double loc[3], std::vector<double>& src)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;
// push function object on stack
    lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(*L, loc[i]);
    lua_pushnumber(*L, tm);
// push field components on stack
    for (unsigned i=0; i<this->getNumInput(); ++i)
      lua_pushnumber(*L, this->getData(i));

// call function
    if (lua_pcall(*L, 4+this->getNumInput(), src.size(), 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("FieldFunctionSource::getSource: ");
      lce << "Problem evaluating function supplied as 'source' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-src.size(); i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("FieldFunctionSource::getSource: Return value not a number");
      src[src.size()+i] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
  }
}
