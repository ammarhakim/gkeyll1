/**
 * @file	LcNodalFiniteElementIfc.cpp
 *
 * @brief       Interface class for a reference nodal finite element.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcGridIfc.h>
#include <LcNodalFiniteElementIfc.h>

namespace Lucee
{
// set module name
  template <> const char *NodalFiniteElementIfc<1>::id = "NodalFiniteElement1D";
  template <> const char *NodalFiniteElementIfc<2>::id = "NodalFiniteElement2D";
  template <> const char *NodalFiniteElementIfc<3>::id = "NodalFiniteElement3D";
  template <> const char *NodalFiniteElementIfc<4>::id = "NodalFiniteElement4D";
  template <> const char *NodalFiniteElementIfc<5>::id = "NodalFiniteElement5D";

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// read in grid
    if (tbl.hasObject<Lucee::GridIfc>("onGrid"))
      grid = &tbl.getObjectAsBase<Lucee::GridIfc>("onGrid");
    else
      throw Lucee::Except("NodalFiniteElementIfc::readInput: Must specify grid using 'onGrid'");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::setIndex(int i) const
  {
    currIdx[0] = i;
  }
  
  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::setIndex(int i, int j) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::setIndex(int i, int j, int k) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
    currIdx[2] = k;
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::setIndex(const int idx[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = idx[i];
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getExclusiveNodeIndices(std::vector<int>& ndIds)
  {
    throw Lucee::Except("NodalFiniteElementIfc::getExclusiveNodeIndices: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  NodalFiniteElementIfc<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNumSurfLowerNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  unsigned
  NodalFiniteElementIfc<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNumSurfUpperNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  unsigned
  NodalFiniteElementIfc<NDIM>::getNumGlobalNodes() const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNumGlobalNodes: Not implemented!");
    return 0;
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfLowerLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfLowerLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfUpperLocalToGlobal(unsigned dim,
    std::vector<int>& lgMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfUpperLocalToGlobal: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfLowerNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfLowerNodeNums: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfUpperNodeNums(unsigned dir,
    std::vector<int>& nodeNum) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfUpperNodeNums: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNodalCoordinates: Not implemented!");    
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getWeights(std::vector<double>& w)
  {
    throw Lucee::Except("NodalFiniteElementIfc::getWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfUpperWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfLowerWeights: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getLowerFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getLowerFaceMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getUpperFaceMassMatrix(unsigned dir,
    Lucee::Matrix<double>& NjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getUpperFaceMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getStiffnessMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getGradStiffnessMatrix(
    unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getGradStiffnessMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  NodalFiniteElementIfc<NDIM>::getNumGaussNodes() const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNumGaussNodes: Not implemented!");
  }

  template <unsigned NDIM>
  unsigned
  NodalFiniteElementIfc<NDIM>::getNumSurfGaussNodes() const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getNumSurfGaussNodes: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfLowerGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getSurfUpperGaussQuadData: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getMomentMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat,
    std::vector<Lucee::Matrix<double> >& lowerMat, std::vector<Lucee::Matrix<double> >& upperMat) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getDiffMatrices: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getHyperDiffusionMatrices(std::vector<Lucee::Matrix<double> >& iMat,
    std::vector<Lucee::Matrix<double> >& lowerMat, std::vector<Lucee::Matrix<double> >& upperMat) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getHyperDiffMatrices: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getLowerReflectingBcMapping(unsigned dir,
    std::vector<unsigned>& nodeMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getLowerReflectingBcMapping: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getUpperReflectingBcMapping(unsigned dir,
    std::vector<unsigned>& nodeMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getLowerReflectingBcMapping: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getLowerFaceToInteriorMapping(unsigned dir,
    Lucee::Matrix<double>& faceToIntMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getLowerFaceToInteriorMapping: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getUpperFaceToInteriorMapping(unsigned dir,
    Lucee::Matrix<double>& faceToIntMap) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getUpperFaceToInteriorMapping: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
    throw Lucee::Except("NodalFiniteElementIfc::extractFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld,
    double *data)
  {
    throw Lucee::Except("NodalFiniteElementIfc::copyAllDataFromField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::copyAllDataToField(const double *data, 
    Lucee::Field<NDIM, double>& fld)
  {
    throw Lucee::Except("NodalFiniteElementIfc::copyAllDataToField: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::evalBasis(double xc[NDIM], std::vector<double>& vals) const;
  {
    throw Lucee::Except("NodalFiniteElementIfc::evalBasis: Not implemented!");
  }

  template <unsigned NDIM>
  NodalFiniteElementIfc<NDIM>::NodalFiniteElementIfc(unsigned numNodes)
    : numNodes(numNodes)
  {
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("numExclusiveNodes", luaNumExclusiveNodes);
    lfm.appendFunc("numNodes", luaNumNodes);
  }

  template <unsigned NDIM>
  int
  NodalFiniteElementIfc<NDIM>::luaNumExclusiveNodes(lua_State *L)
  {
    NodalFiniteElementIfc<NDIM> *f
      = Lucee::PointerHolder<NodalFiniteElementIfc<NDIM> >::getObj(L);

    std::vector<int> ids;
    f->getExclusiveNodeIndices(ids);
    lua_pushnumber(L, ids.size());
    return 1;
  }

  template <unsigned NDIM>
  int
  NodalFiniteElementIfc<NDIM>::luaNumNodes(lua_State *L)
  {
    NodalFiniteElementIfc<NDIM> *f
      = Lucee::PointerHolder<NodalFiniteElementIfc<NDIM> >::getObj(L);
    std::vector<unsigned> ids;
    lua_pushnumber(L, f->getNumNodes());
    return 1;
  }

// instantiations
  template class NodalFiniteElementIfc<1>;
  template class NodalFiniteElementIfc<2>;
  template class NodalFiniteElementIfc<3>;
  template class NodalFiniteElementIfc<4>;
  template class NodalFiniteElementIfc<5>;
}
