/**
 * @file	LcUnstructGridElems.cpp
 *
 * @brief File to define basic elements in unstructured meshes.
 *
 * @version	$Id$
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcUnstructGridElems.h>

namespace Lucee
{
  template <typename REAL>
  void
  VertexElem<REAL>::fillWithCoordinates(REAL xv[3]) const
  {
    for (unsigned i=0; i<3; ++i)
      xv[i] = vcoords[currLoc+i];
  }
  
  template <typename REAL>
  VertexElem<REAL>::VertexElem(const std::vector<REAL>& vc)
    : vcoords(vc), currLoc(0)
  {
  }

  template <typename REAL>
  void
  VertexElem<REAL>::incr() const
  {
    currLoc += 3; // 3 as (x,y,z) coordinates are stored
  }

// instantiations
  template class VertexElem<float>;
  template class VertexElem<double>;
}
