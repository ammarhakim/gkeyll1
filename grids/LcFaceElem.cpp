/**
 * @file	LcFaceElem.cpp
 *
 * @brief       Face elements in unstructured grids.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFaceElem.h>

namespace Lucee
{
  template <typename REAL>
  void
  FaceElem<REAL>::fillWithCoordinates(REAL xv[3]) const
  {
  }

  template <typename REAL>
  REAL
  FaceElem<REAL>::getMeasure() const
  {
    return 0;
  }

  template <typename REAL>
  void
  FaceElem<REAL>::fillWithNormal(REAL norm[3]) const
  {
  }

  template <typename REAL>
  void
  FaceElem<REAL>::fillWithTangent1(REAL tng[3]) const
  {
  }

  template <typename REAL>
  void
  FaceElem<REAL>::fillWithTangent2(REAL tng[3]) const
  {
  }

// instantiations
  template class FaceElem<float>;
  template class FaceElem<double>;
}
