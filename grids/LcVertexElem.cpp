/**
 * @file	LcVertexElem.cpp
 *
 * @brief       Vertex element in unstructured grid.
 */

// lucee includes
#include <LcVertexElem.h>

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

// instantiations
  template class VertexElem<float>;
  template class VertexElem<double>;
}
