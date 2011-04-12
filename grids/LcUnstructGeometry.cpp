/**
 * @file	LcUnstructGeometry.cpp
 *
 * @brief	Class holding geometry of unstructured grids.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructGeometry.h>

namespace Lucee
{
  template <unsigned NDIM, typename REAL>
  UnstructGeometry<NDIM, REAL>::UnstructGeometry()
  {
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumVertices(unsigned nv)
  {
    vcoords.clear();
    vcoords.resize(NDIM*nv);
  }

  template <unsigned NDIM, typename REAL>
  void
  UnstructGeometry<NDIM, REAL>::setNumCells(unsigned nc)
  {
    cellCentroid.clear();
    cellCentroid.resize(NDIM*nc);
    cellVolume.clear();
    cellVolume.resize(nc);
  }

// instantiations
  template class Lucee::UnstructGeometry<1, float>;
  template class Lucee::UnstructGeometry<2, float>;
  template class Lucee::UnstructGeometry<3, float>;

  template class Lucee::UnstructGeometry<1, double>;
  template class Lucee::UnstructGeometry<2, double>;
  template class Lucee::UnstructGeometry<3, double>;
}
