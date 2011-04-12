/**
 * @file	LcEdgeElem.cpp
 *
 * @brief       Edge elements in unstructured grids.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEdgeElem.h>

namespace Lucee
{
  template <typename REAL>
  void
  EdgeElem<REAL>::fillWithCoordinates(REAL xv[3]) const
  {
  }

  template <typename REAL>
  REAL
  EdgeElem<REAL>::getMeasure() const
  {
    return 0;
  }

  template <typename REAL>
  void
  EdgeElem<REAL>::fillWithNormal(REAL norm[3]) const
  {
  }

  template <typename REAL>
  void
  EdgeElem<REAL>::fillWithTangent1(REAL tng[3]) const
  {
  }

  template <typename REAL>
  void
  EdgeElem<REAL>::fillWithTangent2(REAL tng[3]) const
  {
  }

// instantiations
  template class EdgeElem<float>;
  template class EdgeElem<double>;
}
