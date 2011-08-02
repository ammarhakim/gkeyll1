/**
 * @file	LcFaceElem.cpp
 *
 * @brief       Face elements in unstructured grids.
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
    for (unsigned i=0; i<3; ++i)
      xv[i] = fc[3*curr+i];
  }

  template <typename REAL>
  REAL
  FaceElem<REAL>::getMeasure() const
  {
    return fa[curr];
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

  
  template <typename REAL>
  FaceElem<REAL>::FaceElem(const std::vector<REAL>& fc, const std::vector<REAL>& fa)
    : fc(fc), fa(fa), curr(0)
  {
  }

// instantiations
  template class FaceElem<float>;
  template class FaceElem<double>;
}
