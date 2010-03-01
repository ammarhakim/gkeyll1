/**
 * @file	LcField.cpp
 *
 * @brief	Fields hold multiple values per index location.
 *
 * @version	$Id$
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
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(rgn.inflate(0, nc), init),
      numComponents(nc), rgn(rgn), rgnIdx(rgn.inflate(0, nc))
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
    Array<NDIM+1, T, Lucee::RowMajorIndexer>::operator=(val);
    return *this;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>
  Field<NDIM, T>::getView(const Lucee::Region<NDIM, int>& vrgn)
  {
    Array<NDIM+1, T, Lucee::RowMajorIndexer> subArr
      = this->getSlice(vrgn.inflate(0, numComponents));
    Field<NDIM, T> fld(vrgn, numComponents, subArr);
    fld.rgnIdx = this->rgnIdx;
    return fld;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>
  Field<NDIM, T>::getSubCompView(unsigned sc, unsigned ec)
  {
    Array<NDIM+1, T, Lucee::RowMajorIndexer> subArr
      = this->getSlice(rgn.inflate(sc, ec));
    Field<NDIM, T> fld(rgn, ec-sc, subArr);
    fld.rgnIdx = this->rgnIdx;
    int newLower[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i)
      newLower[i] = rgnIdx.getLower(i);
    newLower[NDIM] = -sc; // returned field's 0th component should be sc
    fld.rgnIdx.resetLower(newLower);
    return fld;
  }

  template <unsigned NDIM, typename T>
  Lucee::FieldPtr<T>
  Field<NDIM, T>::createPtr()
  {
    int start[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i)
      start[i] = rgn.getLower(i);
    start[NDIM] = 0;
    unsigned loc = rgnIdx.getIndex(start);
    return Lucee::FieldPtr<T>(numComponents, &this->getRefToLoc(loc));
  }

  template <unsigned NDIM, typename T>
  Lucee::ConstFieldPtr<T>
  Field<NDIM, T>::createConstPtr() const
  {
    int start[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i)
      start[i] = rgn.getLower(i);
    start[NDIM] = 0;
    unsigned loc = rgnIdx.getIndex(start);
    return Lucee::ConstFieldPtr<T>(numComponents, &this->getConstRefToLoc(loc));
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc,
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>& subArr)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(subArr),
      numComponents(nc), rgn(rgn), rgnIdx(rgn.inflate(0, nc))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
  }

// instantiations
  template class Field<1, int>;
  template class Field<2, int>;
  template class Field<3, int>;
  template class Field<4, int>;
  template class Field<5, int>;
  template class Field<6, int>;
  template class Field<7, int>;

  template class Field<1, float>;
  template class Field<2, float>;
  template class Field<3, float>;
  template class Field<4, float>;
  template class Field<5, float>;
  template class Field<6, float>;
  template class Field<7, float>;

  template class Field<1, double>;
  template class Field<2, double>;
  template class Field<3, double>;
  template class Field<4, double>;
  template class Field<5, double>;
  template class Field<6, double>;
  template class Field<7, double>;
}
