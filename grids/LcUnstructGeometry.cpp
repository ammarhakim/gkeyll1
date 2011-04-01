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
  UnstructGeometry<NDIM, REAL>::UnstructGeometry(unsigned nv)
    : vcoords(NDIM*nv)
  {
  }

// instantiations
  template class Lucee::UnstructGeometry<1, float>;
  template class Lucee::UnstructGeometry<2, float>;
  template class Lucee::UnstructGeometry<3, float>;

  template class Lucee::UnstructGeometry<1, double>;
  template class Lucee::UnstructGeometry<2, double>;
  template class Lucee::UnstructGeometry<3, double>;
}
