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
#include <LcFieldFactory.h>

namespace Lucee
{
// names used in registration system
  template <> const char *Field<1, double>::id = "Field1D";
  template <> const char *Field<2, double>::id = "Field2D";
  template <> const char *Field<3, double>::id = "Field3D";
  template <> const char *Field<4, double>::id = "Field4D";
  template <> const char *Field<5, double>::id = "Field5D";
  template <> const char *Field<6, double>::id = "Field6D";
  template <> const char *Field<7, double>::id = "Field7D";

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field()
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(
      &Lucee::FixedVector<NDIM+1, unsigned>(1)[0], (T)0),
      scIdx(0), numComponents(1), rgn(&Lucee::FixedVector<NDIM, int>(1)[0]),
      rgnIdx(rgn.inflate(0, 1))
  {
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(rgn.inflate(0, nc), init),
      scIdx(0), numComponents(nc), rgn(rgn), rgnIdx(rgn.inflate(0, nc))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, 
    int lg[NDIM], int ug[NDIM], const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(
      rgn.extend(lg, ug).inflate(0, nc), init),
      scIdx(0), numComponents(nc), rgn(rgn), rgnIdx(rgn.extend(lg, ug).inflate(0, nc))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = lg[i];
      upperGhost[i] = ug[i];
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Field<NDIM, T>& fld)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(fld),
      scIdx(fld.scIdx),
      numComponents(fld.numComponents),
      rgn(fld.rgn),
      rgnIdx(fld.rgnIdx)
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = fld.lowerGhost[i];
      upperGhost[i] = fld.upperGhost[i];
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::operator=(const Field<NDIM, T>& fld)
  {
    if (&fld == this)
      return *this;

// call base class assignment operator
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>::operator=(fld);
    scIdx = fld.scIdx;
    numComponents = fld.numComponents;
    rgn = fld.rgn;
    rgnIdx = fld.rgnIdx;
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = fld.lowerGhost[i];
      upperGhost[i] = fld.upperGhost[i];
    }
    return *this;
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
  void
  Field<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::FieldFactory<NDIM, T> ff;
    ff.readInput(tbl);
    *this = *ff.create();
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>
  Field<NDIM, T>::getView(const Lucee::Region<NDIM, int>& vrgn)
  {
    Array<NDIM+1, T, Lucee::RowMajorIndexer> subArr
      = this->getSlice(
        vrgn.extend(lowerGhost, upperGhost).inflate(scIdx, numComponents+scIdx));
    Field<NDIM, T> fld(vrgn, scIdx, scIdx+numComponents, lowerGhost, upperGhost, subArr);
    fld.rgnIdx = this->rgnIdx;
    return fld;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>
  Field<NDIM, T>::getSubCompView(unsigned sc, unsigned ec)
  {
    Array<NDIM+1, T, Lucee::RowMajorIndexer> subArr
      = this->getSlice(
        rgn.extend(lowerGhost, upperGhost).inflate(sc, ec));

    int newLower[NDIM+1];
    for (unsigned i=0; i<NDIM; ++i)
      newLower[i] = rgnIdx.getLower(i);
    newLower[NDIM] = -sc; // returned field's 0th component should be sc

    subArr.resetLower(newLower);
    Field<NDIM, T> fld(rgn, sc, ec, lowerGhost, upperGhost, subArr);

    fld.rgnIdx = this->rgnIdx;
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
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned sc, unsigned ec,
    int lg[NDIM], int ug[NDIM], Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>& subArr)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(subArr),
      scIdx(sc), numComponents(ec-sc), rgn(rgn), rgnIdx(rgn.inflate(sc, ec))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = lg[i];
      upperGhost[i] = ug[i];
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
