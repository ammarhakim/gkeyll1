/**
 * @file	LcNodalDgFunctionBoundaryCondition.cpp
 *
 * @brief	Class for applying copy boundary conditions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcLuaState.h>
#include <LcNodalDgFunctionBoundaryCondition.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set costructor name
  template <> const char *NodalDgFunctionBoundaryCondition<1>::id = "NodalDgFunction1D";
  template <> const char *NodalDgFunctionBoundaryCondition<2>::id = "NodalDgFunction2D";
  template <> const char *NodalDgFunctionBoundaryCondition<3>::id = "NodalDgFunction3D";

  template <unsigned NDIM>
  void
  NodalDgFunctionBoundaryCondition<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    BoundaryCondition::readInput(tbl);
// get reference to function
    fnRef = tbl.getFunctionRef("bc");
    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("NodalDisContHyperUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  NodalDgFunctionBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
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
    if (lua_pcall(*L, 4, this->numComponents(), 0) != 0)
    {
      std::string err(lua_tostring(*L, -1));
      lua_pop(*L, 1);
      Lucee::Except lce("NodalDgFunctionBoundaryCondition::applyBc: ");
      lce << "Problem evaluating function supplied as 'bc' ";
      lce << std::endl << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-this->numComponents(); i<0; ++i)
    {
      if (!lua_isnumber(*L, i))
        throw Lucee::Except("NodalDgFunctionBoundaryCondition::applyBc: Return value not a number");
      qbc[this->component(this->numComponents()+i)] = lua_tonumber(*L, i);
    }
    lua_pop(*L, 1);
  }

// instatiations
  template class NodalDgFunctionBoundaryCondition<1>;
  template class NodalDgFunctionBoundaryCondition<2>;
  template class NodalDgFunctionBoundaryCondition<3>;
}
