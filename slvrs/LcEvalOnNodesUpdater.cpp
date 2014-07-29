/**
 * @file	LcEvalOnNodesUpdater.cpp
 *
 * @brief	Set region based on predicate.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEvalOnNodesUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{

  template <> const char *EvalOnNodesUpdater<1>::id = "EvalOnNodes1D";
  template <> const char *EvalOnNodesUpdater<2>::id = "EvalOnNodes2D";
  template <> const char *EvalOnNodesUpdater<3>::id = "EvalOnNodes3D";
  template <> const char *EvalOnNodesUpdater<4>::id = "EvalOnNodes4D";
  template <> const char *EvalOnNodesUpdater<5>::id = "EvalOnNodes5D";

  template <unsigned NDIM>
  EvalOnNodesUpdater<NDIM>::EvalOnNodesUpdater()
    : fnRef(-1), sharedNodes(false)
  {
  }

  template <unsigned NDIM>
  void
  EvalOnNodesUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
// get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    sharedNodes = false;
// check if there are shared nodes
    if (tbl.hasBool("shareCommonNodes"))
      sharedNodes = tbl.getBool("shareCommonNodes");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("EvalOnNodesUpdater::readInput: Must specify element to use using 'basis'");
  }

  template <unsigned NDIM>
  void
  EvalOnNodesUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  EvalOnNodesUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);
    q = 0.0;

    std::vector<int> ndIds;
    if (sharedNodes)
      nodalBasis->getExclusiveNodeIndices(ndIds);
    else
    { // create "unit" mapping
      ndIds.resize(nodalBasis->getNumNodes());
      for (unsigned i = 0; i<ndIds.size(); ++i)
        ndIds[i] = i;
    }

    unsigned numNodes = ndIds.size();
    unsigned nc = q.getNumComponents()/numNodes;
    std::vector<double> res(nc);
    int idx[NDIM];

    int numCoordinates = NC;
    Lucee::Matrix<double> nodeCoords(nodalBasis->getNumNodes(), numCoordinates);

    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    Lucee::FieldPtr<double> ptr = q.createPtr();
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);

      nodalBasis->setIndex(idx);
      nodalBasis->getNodalCoordinates(nodeCoords);

      for (unsigned n=0; n<numNodes; ++n)
      {
        evaluateFunction(*L, t, nodeCoords, ndIds[n], res);
        for (unsigned k=0; k<nc; ++k)
          ptr[nc*n+k] = res[k];
      }
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  EvalOnNodesUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  EvalOnNodesUpdater<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    const Lucee::Matrix<double> nc, unsigned nn, std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    for (unsigned i=0; i<NC; ++i)
      lua_pushnumber(L, nc(nn,i));
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, NC+1, res.size(), 0) != 0)
    {
      Lucee::Except lce("EvalOnNodesUpdater::evaluateFunction: ");
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
        throw Lucee::Except("EvalOnNodesUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class EvalOnNodesUpdater<1>;
  template class EvalOnNodesUpdater<2>;
  template class EvalOnNodesUpdater<3>;
  template class EvalOnNodesUpdater<4>;
  template class EvalOnNodesUpdater<5>;
}
