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
  NodalFiniteElementIfc<NDIM>::getMassMatrix(Lucee::Matrix<double> NjNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getMassMatrix: Not implemented!");
  }

  template <unsigned NDIM>
  void
  NodalFiniteElementIfc<NDIM>::getStiffnessMatrix(Lucee::Matrix<double> DNjDNk) const
  {
    throw Lucee::Except("NodalFiniteElementIfc::getStiffnessMatrix: Not implemented!");
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
  NodalFiniteElementIfc<NDIM>::NodalFiniteElementIfc(unsigned numNodes)
    : numNodes(numNodes)
  {
  }

// instantiations
  template class NodalFiniteElementIfc<1>;
  template class NodalFiniteElementIfc<2>;
  template class NodalFiniteElementIfc<3>;
}
