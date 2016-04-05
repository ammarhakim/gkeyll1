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
      throw Lucee::Except("NodalDgFunctionBoundaryCondition::readInput: Must specify element to use using 'basis'");
// matrix to store coordinates of nodes
    nodeCoords = Lucee::Matrix<double>(nodalBasis->getNumNodes(), 3);
// create unit mapping
    ndIds.resize(nodalBasis->getNumNodes());
    for (unsigned i=0; i<ndIds.size(); ++i) ndIds[i] = i;
  }

  template <unsigned NDIM>
  void
  NodalDgFunctionBoundaryCondition<NDIM>::applyBc(double tm, const double loc[3], const int *idx,
    const Lucee::RectCoordSys& c, const Lucee::ConstFieldPtr<double>& qin1,
    const Lucee::ConstFieldPtr<double>& qin, Lucee::FieldPtr<double>& qbc)
  {
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    unsigned numNodes = nodalBasis->getNumNodes();
    unsigned nc = qbc.getNumComponents()/numNodes;
    std::vector<double> res(nc);

// get nodal coordinates
    nodalBasis->setIndex(idx);
    nodalBasis->getNodalCoordinates(nodeCoords);

    for (unsigned n=0; n<numNodes; ++n)
    {
      evaluateFunction(*L, tm, nodeCoords, ndIds[n], res);
      for (unsigned k=0; k<nc; ++k)
        qbc[nc*n+k] = res[k];
    }
  }

  template <unsigned NDIM>
  void
  NodalDgFunctionBoundaryCondition<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    const Lucee::Matrix<double>& nc, unsigned nn, std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<3; ++i)
      lua_pushnumber(L, nc(nn,i));
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, 4, res.size(), 0) != 0)
    {
      Lucee::Except lce("NodalDgFunctionBoundaryCondition::evaluateFunction: ");
      lce << "Problem evaluating function supplied as 'evaluate' "
          << std::endl;
      std::string err(lua_tostring(L, -1));
      lua_pop(L, 1);
      lce << "[" << err << "]";
      throw lce;
    }
// fetch results
    for (int i=-res.size(); i<0; ++i)
    {
      if (!lua_isnumber(L, i))
        throw Lucee::Except("NodalDgFunctionBoundaryCondition::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instatiations
  template class NodalDgFunctionBoundaryCondition<1>;
  template class NodalDgFunctionBoundaryCondition<2>;
  template class NodalDgFunctionBoundaryCondition<3>;
}
