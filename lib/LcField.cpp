/**
 * @file	LcField.cpp
 *
 * @brief	Fields hold multiple values per index location.
 *
 * @version	$Id: LcField.cpp 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>

namespace Lucee
{
  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >(rgn.inflate(0, nc), init),
      numComponents(nc), rgn(rgn), rgnIdx(rgn)
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::operator=(const T& val)
  {
// simply call base class assignment operator    
    Array<NDIM+1, T, Lucee::RowMajorIndexer<NDIM+1> >::operator=(val);
    return *this;
  }

  template <unsigned NDIM, typename T>
  Lucee::FieldPtr<T>
  Field<NDIM, T>::createPtr()
  {
    int start[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
      start[i] = rgn.getLower(i);
    unsigned loc = rgnIdx.getGenIndex(start);
    return Lucee::FieldPtr<T>(numComponents, &this->getRefToLoc(loc));
  }

// instantiations
  template class Lucee::Field<1, int>;
  template class Lucee::Field<1, float>;
  template class Lucee::Field<1, double>;

  template class Lucee::Field<2, int>;
  template class Lucee::Field<2, float>;
  template class Lucee::Field<2, double>;

  template class Lucee::Field<3, int>;
  template class Lucee::Field<3, float>;
  template class Lucee::Field<3, double>;
}
