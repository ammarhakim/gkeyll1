/**
 * @file	LcLagrangeTensorElement.cpp
 *
 * @brief       Lagrange tensor-product element.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLagrangeTensorElement.h>

namespace Lucee
{
  template <> const char *LagrangeTensorElement<1>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<2>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<3>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<4>::id = "LagrangeTensor";
  template <> const char *LagrangeTensorElement<5>::id = "LagrangeTensor";

  template <unsigned NDIM>
  LagrangeTensorElement<NDIM>::LagrangeTensorElement()
    : Lucee::NodalFiniteElementIfc<NDIM>(1)
  {
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// call base class method
    Lucee::NodalFiniteElementIfc<NDIM>::readInput(tbl);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfLowerNodes(unsigned dir) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumSurfLowerNodes(dir);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfUpperNodes(unsigned dir) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumSurfUpperNodes(dir);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumGlobalNodes() const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumGlobalNodes();
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getExclusiveNodeIndices(std::vector<unsigned>& ndIds)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getExclusiveNodeIndices(ndIds);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getLocalToGlobal(std::vector<int>& lgMap) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getLocalToGlobal(lgMap);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerLocalToGlobal(dir, lgMap);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperLocalToGlobal(unsigned dir, std::vector<int>& lgMap) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfUpperLocalToGlobal(dir, lgMap);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerNodeNums(dir, nodeNum);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperNodeNums(unsigned dir, std::vector<int>& nodeNum) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfUpperNodeNums(dir, nodeNum);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getNodalCoordinates(Lucee::Matrix<double>& nodeCoords)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNodalCoordinates(nodeCoords);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getWeights(std::vector<double>& w)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getWeights(w);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperWeights(unsigned dir, std::vector<double>& w)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfUpperWeights(dir, w);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerWeights(unsigned dir, std::vector<double>& w)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerWeights(dir, w);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getMassMatrix(Lucee::Matrix<double>& NjNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getMassMatrix(NjNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getLowerFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& NjNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getLowerFaceMassMatrix(dir, NjNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getUpperFaceMassMatrix(unsigned dir, Lucee::Matrix<double>& NjNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getUpperFaceMassMatrix(dir, NjNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getStiffnessMatrix(Lucee::Matrix<double>& DNjDNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getStiffnessMatrix(DNjDNk);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getGradStiffnessMatrix(unsigned dir, Lucee::Matrix<double>& DNjNk) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getGradStiffnessMatrix(dir, DNjNk);
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumGaussNodes() const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumGaussNodes();
  }

  template <unsigned NDIM>
  unsigned
  LagrangeTensorElement<NDIM>::getNumSurfGaussNodes() const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getNumSurfGaussNodes();
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getGaussQuadData(Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getGaussQuadData(interpMat, ordinates, weights);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfLowerGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getSurfLowerGaussQuadData(dir, interpMat, ordinates, weights);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getSurfUpperGaussQuadData(unsigned dir, Lucee::Matrix<double>& interpMat,
    Lucee::Matrix<double>& ordinates, std::vector<double>& weights) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>:: getSurfUpperGaussQuadData(dir, interpMat, ordinates, weights);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::getMomentMatrix(unsigned p, Lucee::Matrix<double>& momMat) const
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::getMomentMatrix(p, momMat);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::extractFromField(const Lucee::Field<NDIM, double>& fld,
    std::vector<double>& data)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::extractFromField(fld, data);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataFromField(const Lucee::Field<NDIM, double>& fld, double *data)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::copyAllDataFromField(fld, data);
  }

  template <unsigned NDIM>
  void
  LagrangeTensorElement<NDIM>::copyAllDataToField(const double *data, Lucee::Field<NDIM, double>& fld)
  {
    return Lucee::NodalFiniteElementIfc<NDIM>::copyAllDataToField(data, fld);
  }

// instantiations
  template class LagrangeTensorElement<1>;
  template class LagrangeTensorElement<2>;
  template class LagrangeTensorElement<3>;
  template class LagrangeTensorElement<4>;
  template class LagrangeTensorElement<5>;
}
