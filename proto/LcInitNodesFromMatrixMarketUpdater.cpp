/**
 * @file	LcInitNodesFromMatrixMarketUpdater.cpp
 *
 * @brief	Given an input file (.mtx), reads data into a field
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcInitNodesFromMatrixMarketUpdater.h>
#include <LcGlobals.h>

// loki includes
#include <loki/Singleton.h>

// std includes
#include <vector>

// eigen include to write to file
#include <unsupported/Eigen/SparseExtra>

namespace Lucee
{

  template <> const char *InitNodesFromMatrixMarketUpdater<1>::id = "InitNodesFromMatrixMarket1D";
  template <> const char *InitNodesFromMatrixMarketUpdater<2>::id = "InitNodesFromMatrixMarket2D";
  template <> const char *InitNodesFromMatrixMarketUpdater<3>::id = "InitNodesFromMatrixMarket3D";
  template <> const char *InitNodesFromMatrixMarketUpdater<4>::id = "InitNodesFromMatrixMarket4D";
  template <> const char *InitNodesFromMatrixMarketUpdater<5>::id = "InitNodesFromMatrixMarket5D";

  template <unsigned NDIM>
  InitNodesFromMatrixMarketUpdater<NDIM>::InitNodesFromMatrixMarketUpdater()
    : UpdaterIfc()
  {
  }

  template <unsigned NDIM>
  void
  InitNodesFromMatrixMarketUpdater<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
    UpdaterIfc::readInput(tbl);

    if (tbl.hasObject<Lucee::NodalFiniteElementIfc<NDIM> >("basis"))
      nodalBasis = &tbl.getObjectAsBase<Lucee::NodalFiniteElementIfc<NDIM> >("basis");
    else
      throw Lucee::Except("InitNodesFromMatrixMarketUpdater::readInput: Must specify element to use using 'basis'");

    if (tbl.hasNumber("totalNodes"))
      totalNodes = tbl.getNumber("totalNodes");
    else
      throw Lucee::Except("InitNodesFromMatrixMarketUpdater::readInput: Must specify totalNodes");

    if (tbl.hasString("filename"))
      filename = tbl.getString("filename");
    else
      throw Lucee::Except("InitNodesFromMatrixMarketUpdater::readInput: Must specify filename");
      
  }

  template <unsigned NDIM>
  void
  InitNodesFromMatrixMarketUpdater<NDIM>::initialize()
  {
    UpdaterIfc::initialize();

    // Read data
    Eigen::loadMarket(initialValueMatrix, filename);
  }

  template <unsigned NDIM>
  Lucee::UpdaterStatus
  InitNodesFromMatrixMarketUpdater<NDIM>::update(double t)
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
    //Lucee::RowMajorSequencer<NDIM> seq(localRgn);
    Lucee::Region<NDIM, int> localExtRgn = q.getExtRegion();
    Lucee::RowMajorSequencer<NDIM> seq(localExtRgn);
    Lucee::RowMajorIndexer<NDIM> volIdxr(localExtRgn);
    // Figure out what cell index to set to zero
    //evaluateFunction(*L, t, res);
    // This will be a number from 0 to NodesPerCell*TotalCells
    //int targetIndex = (int) res[0];
    // This is the actual index of the node within a cell
    //int targetNode = targetIndex % nlocal;
    
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      q.setPtr(ptr, idx);
      nodalBasis->setIndex(idx);

      int cellIndex = volIdxr.getIndex(idx);

      for (int nodeIndex = 0; nodeIndex < nlocal; nodeIndex++)
        ptr[nodeIndex] = initialValueMatrix.coeffRef(cellIndex*nlocal+nodeIndex,0);
    }

    return Lucee::UpdaterStatus();
  }

  template <unsigned NDIM>
  void
  InitNodesFromMatrixMarketUpdater<NDIM>::declareTypes()
  {
    this->appendOutVarType(typeid(Lucee::Field<NDIM, double>));
  }

// instantiations
  template class InitNodesFromMatrixMarketUpdater<1>;
  template class InitNodesFromMatrixMarketUpdater<2>;
  template class InitNodesFromMatrixMarketUpdater<3>;
  template class InitNodesFromMatrixMarketUpdater<4>;
  template class InitNodesFromMatrixMarketUpdater<5>;
}
