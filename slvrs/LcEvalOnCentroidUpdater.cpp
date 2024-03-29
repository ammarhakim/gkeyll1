/**
 * @file	LcEvalOnCentroidUpdater.cpp
 *
 * @brief	Set region based on predicate.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEvalOnCentroidUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{

  template <> const char *EvalOnCentroidUpdater<1>::id = "EvalOnCentroid1D";
  template <> const char *EvalOnCentroidUpdater<2>::id = "EvalOnCentroid2D";
  template <> const char *EvalOnCentroidUpdater<3>::id = "EvalOnCentroid3D";

  template <unsigned NDIM>
  EvalOnCentroidUpdater<NDIM>::EvalOnCentroidUpdater()
    : fnRef(-1), sharedCentroid(false)
  {
  }

  template <unsigned NDIM>
  void
  EvalOnCentroidUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
// get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("EvalOnNodesUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  EvalOnCentroidUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  EvalOnCentroidUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);
    q = 0.0;

    unsigned numNodes = nodalBasis->getNumNodes();
    unsigned nc = q.getNumComponents()/numNodes;
    std::vector<double> res(nc);
    int idx[NDIM];

    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    double xc[3];
    Lucee::FieldPtr<double> ptr = q.createPtr();
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);
      grid.setIndex(idx);
      grid.getCentroid(xc);

      nodalBasis->setIndex(idx);

      for (unsigned n=0; n<numNodes; ++n)
      {
        evaluateFunction(*L, t, xc, 1, res);
        for (unsigned k=0; k<nc; ++k)
          ptr[nc*n+k] = res[k];
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  EvalOnCentroidUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  EvalOnCentroidUpdater<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    double nc[3], unsigned nn, std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<NC; ++i)
      lua_pushnumber(L, nc[i]);
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, NC+1, res.size(), 0) != 0)
    {
      Lucee::Except lce("EvalOnCentroidUpdater::evaluateFunction: ");
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
        throw Lucee::Except("EvalOnCentroidUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class EvalOnCentroidUpdater<1>;
  template class EvalOnCentroidUpdater<2>;
  template class EvalOnCentroidUpdater<3>;
}
