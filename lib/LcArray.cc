/**
 * @file	LcArray.cc
 *
 * @brief	Serial array class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcArray.h>

namespace Lucee
{
  template <unsigned NDIM, typename T>
  Array<NDIM, T>::Array(unsigned shp[NDIM], const T& init) 
    : Lucee::Object("Lucee::Array")
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = shp[i];
      start[i] = 0;
    }
  }

  template <unsigned NDIM, typename T>
  Array<NDIM, T>::Array(unsigned shp[NDIM], int sta[NDIM], const T& init)
  : Lucee::Object("Lucee::Array")
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      shape[i] = shp[i];
      start[i] = sta[i];
    }
  }

  template <unsigned NDIM, typename T>
  void 
  Array<NDIM, T>::getShape(unsigned shp[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      shp[i] = shape[i];
  }

  template <unsigned NDIM, typename T>
  unsigned 
  Array<NDIM, T>::getShape(unsigned dir) const
  {
    return shape[dir];
  }

  template <unsigned NDIM, typename T>
  int
  Array<NDIM, T>::getStart(unsigned dir) const
  {
    return start[dir];
  }

  template <unsigned NDIM, typename T>
  int
  Array<NDIM, T>::getEnd(unsigned dir) const
  {
    return start[dir]+shape[dir];
  }

// instantiations
  template class Array<1, int>;
  template class Array<2, int>;
  template class Array<3, int>;
  template class Array<4, int>;

  template class Array<1, float>;
  template class Array<2, float>;
  template class Array<3, float>;
  template class Array<4, float>;

  template class Array<1, double>;
  template class Array<2, double>;
  template class Array<3, double>;
  template class Array<4, double>;
}
