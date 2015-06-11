/**
 * @file	LcSetSingleNodeToOneUpdater.cpp
 *
 * @brief	Set region based on predicate.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSetSingleNodeToOneUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

namespace Lucee
{

  template <> const char *SetSingleNodeToOneUpdater<1>::id = "SetSingleNodeToOne1D";
  template <> const char *SetSingleNodeToOneUpdater<2>::id = "SetSingleNodeToOne2D";
  template <> const char *SetSingleNodeToOneUpdater<3>::id = "SetSingleNodeToOne3D";
  template <> const char *SetSingleNodeToOneUpdater<4>::id = "SetSingleNodeToOne4D";
  template <> const char *SetSingleNodeToOneUpdater<5>::id = "SetSingleNodeToOne5D";

  template <unsigned NDIM>
  SetSingleNodeToOneUpdater<NDIM>::SetSingleNodeToOneUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  SetSingleNodeToOneUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);
    // get function to evaluate
    fnRef = tbl.getFunctionRef("evaluate");

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("SetSingleNodeToOneUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("scaleFactor"))
      scaleFactor = tbl.getNumber("scaleFactor");
    else scaleFactor = 1.0;
  }

  template <unsigned NDIM>
  void
  SetSingleNodeToOneUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  SetSingleNodeToOneUpdater<NDIM>::update(double t)
  {
    const Lucee::StructuredGridBase<NDIM>& grid 
      = this->getGrid<Lucee::StructuredGridBase<NDIM> >();
    Lucee::Field<NDIM, double>& q = this->getOut<Lucee::Field<NDIM, double> >(0);
    q = 0.0;

    std::vector<double> res(1);
    int idx[NDIM];
    int nlocal = nodalBasis->getNumNodes();
    
    Lucee::LuaState *L = Loki::SingletonHolder<Lucee::Globals>::Instance().L;

    Lucee::FieldPtr<double> ptr = q.createPtr();
    Lucee::Region<NDIM, int> localRgn = grid.getLocalRegion();
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();

    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    
    Lucee::RowMajorIndexer<NDIM> volIdxr(localExtRgn);
    // Figure out what cell index to set to zero
    evaluateFunction(*L, t, res);
    // This will be a number from 0 to NodesPerCell*TotalCells
    int targetIndex = (int) res[0];
    // This is the actual index of the node within a cell
    int targetNode = targetIndex % nlocal;

    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);
      nodalBasis->setIndex(idx);

      int cellIndex = volIdxr.getIndex(idx);

      // If the target node is in this cell, then set its value to 1.0
      if (targetIndex < (cellIndex+1)*nlocal && targetIndex >= cellIndex*nlocal)
        ptr[targetNode] = scaleFactor;
      else
        ptr[targetNode] = 0.0;
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  SetSingleNodeToOneUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

  template <unsigned NDIM>
  void
  SetSingleNodeToOneUpdater<NDIM>::evaluateFunction(Lucee::LuaState& L, double tm,
    std::vector<double>& res)
  {
// push function object on stack
    lua_rawgeti(L, LUA_REGISTRYINDEX, fnRef);
// push variables on stack
    lua_pushnumber(L, tm);
// call function
    if (lua_pcall(L, 1, res.size(), 0) != 0)
    {
      Lucee::Except lce("SetSingleNodeToOneUpdater::evaluateFunction: ");
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
        throw Lucee::Except("SetSingleNodeToOneUpdater::evaluateFunction: Return value not a number");
      res[res.size()+i] = lua_tonumber(L, i);
    }
    lua_pop(L, 1);
  }

// instantiations
  template class SetSingleNodeToOneUpdater<1>;
  template class SetSingleNodeToOneUpdater<2>;
  template class SetSingleNodeToOneUpdater<3>;
  template class SetSingleNodeToOneUpdater<4>;
  template class SetSingleNodeToOneUpdater<5>;
}
