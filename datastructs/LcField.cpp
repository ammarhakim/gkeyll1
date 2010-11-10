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
#include <LcPointerHolder.h>

namespace Lucee
{
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

    Field<NDIM, T> fld(rgn, sc, ec, lowerGhost, upperGhost, subArr);
    fld.resetLowerForIndexer(newLower); // reset the Array indexer
    fld.rgnIdx = this->rgnIdx;
    fld.rgnIdx.resetLower(newLower);

    newLower[NDIM] = 0; // start index should be 0, always
    fld.resetLower(newLower); // reset Array start

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

  template <unsigned NDIM, typename T>
  Lucee::IoNodeType
  Field<NDIM, T>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node, const std::string& nm)
  {
    std::vector<size_t> dataSetSize(NDIM+1), dataSetBeg(NDIM+1), dataSetLen(NDIM+1);
// construct sizes and shapes to write stuff out
    for (unsigned i=0; i<NDIM; ++i)
    {
      unsigned dirShape = getUpper(i)-getLower(i);
      dataSetSize[i] = dirShape;
      dataSetBeg[i] = 0;
      dataSetLen[i] = dirShape;
    }
    dataSetSize[NDIM] = numComponents;
    dataSetBeg[NDIM] = 0;
    dataSetLen[NDIM] = numComponents;

    Lucee::Region<NDIM+1, int> myRgn = rgn.inflate(0, numComponents);
    std::vector<T> buff(myRgn.getVolume());
    Lucee::RowMajorSequencer<NDIM+1> seq(myRgn); // must be row-major for HDF5
// copy data into buffer
    unsigned count = 0;
    while (seq.step())
      buff[count++] = this->operator()(seq.getIndex());
// write it out
    Lucee::IoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    return dn;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::copy(const Field<NDIM, T>& fld)
  {
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>::copy(fld);
    return *this;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::accumulate(double coeff, const Field<NDIM, T>& fld)
  {
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>::accumulate(coeff, fld);
    return *this;
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("clear", luaClear);
    lfm.appendFunc("copy", luaCopy);
    lfm.appendFunc("accumulate", luaAccumulate);
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaClear(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (! lua_isnumber(L, 2))
    {
      Lucee::Except lce("Field::luaClear: Must provide a number to 'clear' method");
      throw lce;
    }
    T num = (T) lua_tonumber(L, 2);
    (*fld) = num;
    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaCopy(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (lua_type(L, 2) != LUA_TUSERDATA)
    {
      Lucee::Except lce("Field::luaCopy: Must provide a field to 'copy' method");
      throw lce;
    }
    Lucee::PointerHolder<Field<NDIM, T> > *fldPtr =
      (Lucee::PointerHolder<Field<NDIM, T> >*) lua_touserdata(L, 2);
    fld->copy(*fldPtr->pointer);

    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaAccumulate(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (! lua_isnumber(L, 2))
    {
      Lucee::Except lce("Field::luaAccumulate: Must provide a number to 'accumulatr' method");
      throw lce;
    }
    T coeff = (T) lua_tonumber(L, 2);
    if (lua_type(L, 3) != LUA_TUSERDATA)
    {
      Lucee::Except lce("Field::luaAccumulate: Must provide a field to 'accumulate' method");
      throw lce;
    }
    Lucee::PointerHolder<Field<NDIM, T> > *fldPtr =
      (Lucee::PointerHolder<Field<NDIM, T> >*) lua_touserdata(L, 3);
    fld->accumulate(coeff, *fldPtr->pointer);

    return 0;
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
