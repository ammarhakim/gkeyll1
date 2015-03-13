/**
 * @file	LcField.cpp
 *
 * @brief	Fields hold multiple values per index location.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcField.h>
#include <LcGlobals.h>
#include <LcMathLib.h>
#include <LcPointerHolder.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field()
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(
      &Lucee::FixedVector<NDIM+1, unsigned>(1)[0], (T)0),
      scIdx(0), numComponents(1), rgn(&Lucee::FixedVector<NDIM, int>(1)[0]),
      rgnIdx(rgn.inflate(0, 1))
  {
    globalRgn = rgn; // for now global region = local region
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& rgn, unsigned nc, const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(rgn.inflate(0, nc), init),
      scIdx(0), numComponents(nc), rgn(rgn), rgnIdx(rgn.inflate(0, nc))
  {
    globalRgn = rgn;  // for now global region = local region
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
    globalRgn = rgn;  // for now global region = local region
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = lg[i];
      upperGhost[i] = ug[i];
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& globalRgn, const Lucee::Region<NDIM, int>& rgn,
    unsigned nc, const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(rgn.inflate(0, nc), init),
      scIdx(0), numComponents(nc), rgn(rgn), globalRgn(globalRgn), rgnIdx(rgn.inflate(0, nc))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& globalRgn, const Lucee::Region<NDIM, int>& rgn,
    unsigned nc, int lg[NDIM], int ug[NDIM], const T& init)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(
      rgn.extend(lg, ug).inflate(0, nc), init),
      scIdx(0), numComponents(nc), rgn(rgn), globalRgn(globalRgn), rgnIdx(rgn.extend(lg, ug).inflate(0, nc))
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
      rgnIdx(fld.rgnIdx),
      globalRgn(fld.globalRgn)
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
    globalRgn = fld.globalRgn;
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
    Field<NDIM, T> fld(globalRgn, vrgn, scIdx, scIdx+numComponents, lowerGhost, upperGhost, subArr);
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

    Field<NDIM, T> fld(globalRgn, rgn, sc, ec, lowerGhost, upperGhost, subArr);
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
  Field<NDIM, T>::Field(const Lucee::Region<NDIM, int>& globalRgn, 
    const Lucee::Region<NDIM, int>& rgn, unsigned sc, unsigned ec,
    int lg[NDIM], int ug[NDIM], Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>& subArr)
    : Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>(subArr),
      scIdx(sc), numComponents(ec-sc), rgn(rgn), globalRgn(globalRgn), rgnIdx(rgn.inflate(sc, ec))
  {
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = lg[i];
      upperGhost[i] = ug[i];
    }
  }

  template <unsigned NDIM, typename T>
  TxIoNodeType
  Field<NDIM, T>::writeToFile(TxIoBase& io, TxIoNodeType& node, const std::string& nm)
  {
    std::vector<size_t> dataSetSize(NDIM+1), dataSetBeg(NDIM+1), dataSetLen(NDIM+1);
// construct sizes and shapes to write stuff out
    for (unsigned i=0; i<NDIM; ++i)
    {
      dataSetSize[i] = globalRgn.getShape(i);
      dataSetBeg[i] = rgn.getLower(i);
      dataSetLen[i] = rgn.getShape(i);
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
    TxIoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    return dn;
  }

  template <unsigned NDIM, typename T>
  TxIoNodeType
  Field<NDIM, T>::writeToFileWithGhost(TxIoBase& io, TxIoNodeType& node, 
    int lg[NDIM], int ug[NDIM], const std::string& nm)
  {
    Lucee::Region<NDIM, int> globalRgnGst = globalRgn.extend(lg, ug);
    Lucee::Region<NDIM, int> rgnGst = rgn.extend(lg, ug);

    std::vector<size_t> dataSetSize(NDIM+1), dataSetBeg(NDIM+1), dataSetLen(NDIM+1);
// construct sizes and shapes to write stuff out
    for (unsigned i=0; i<NDIM; ++i)
    {
      dataSetSize[i] = globalRgnGst.getShape(i);
// HDF5 expects global start index is always zero-based
      dataSetBeg[i] = rgnGst.getLower(i) - globalRgnGst.getLower(i);
      dataSetLen[i] = rgnGst.getShape(i);
    }
    dataSetSize[NDIM] = numComponents;
    dataSetBeg[NDIM] = 0;
    dataSetLen[NDIM] = numComponents;

    Lucee::Region<NDIM+1, int> myRgn = rgnGst.inflate(0, numComponents);
    std::vector<T> buff(myRgn.getVolume());
    Lucee::RowMajorSequencer<NDIM+1> seq(myRgn); // must be row-major for HDF5
// copy data into buffer
    unsigned count = 0;
    while (seq.step())
      buff[count++] = this->operator()(seq.getIndex());
// write it out
    TxIoNodeType dn =
      io.writeDataSet(node, nm, dataSetSize, dataSetBeg, dataSetLen, &buff[0]);

    return dn;
  }

    template <unsigned NDIM, typename T>
    TxIoNodeType
    Field<NDIM, T>::readFromFile(TxIoBase& io, TxIoNodeType& node, const std::string& nm)
    {
      std::vector<size_t> dataSetSize(NDIM+1), dataSetBeg(NDIM+1), dataSetLen(NDIM+1);

// construct sizes and shapes to read stuff in
      for (unsigned i=0; i<NDIM; ++i)
      {
        dataSetSize[i] = globalRgn.getShape(i);
        dataSetBeg[i] = rgn.getLower(i);
        dataSetLen[i] = rgn.getShape(i);
      }
      dataSetSize[NDIM] = numComponents;
      dataSetBeg[NDIM] = 0;
      dataSetLen[NDIM] = numComponents;

      Lucee::Region<NDIM+1, int> myRgn = rgn.inflate(0, numComponents);
      std::vector<T> rdata(myRgn.getVolume());
      TxIoNodeType dn = io.readDataSet(node, nm, dataSetBeg, dataSetLen, &rdata[0]);

      Lucee::RowMajorSequencer<NDIM+1> seq(myRgn); // must be row-major for HDF5
// copy data from buffer
      unsigned count = 0;
      while (seq.step())
        this->operator()(seq.getIndex()) = rdata[count++];

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
  Field<NDIM, T>
  Field<NDIM, T>::duplicate()
  {
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer> arr 
      = Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>::duplicate();
    Lucee::Field<NDIM, T> fld(globalRgn, rgn, 0, numComponents, lowerGhost, upperGhost, arr);
    fld.rgnIdx = this->rgnIdx;
    return fld;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::accumulate(double coeff, const Field<NDIM, T>& fld)
  {
    Lucee::Array<NDIM+1, T, Lucee::RowMajorIndexer>::accumulate(coeff, fld);
    return *this;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::applyPeriodicBc(unsigned dir)
  {
// DOES NOT APPLY PERIODIC BCS TO CORNER CELLS. ALSO NEED TO DO FUNKY STUFF FOR PARALLEL

    int lo[NDIM], up[NDIM];

// create a region to represent lower interior layer of cells
    for (unsigned i=0; i<NDIM; ++i)
    { // whole interior region
      lo[i] = getLower(i);
      up[i] = getUpper(i);
    }
// adjust cells along 'dir' direction so it represents region to be copied
    up[dir] = lo[dir]+upperGhost[dir];
    Lucee::Region<NDIM, int> lowerRgn(lo, up);

    unsigned stride = rgn.getShape(dir); // no of cells in interior in specified direction

    Lucee::ConstFieldPtr<T> ptr = this->createConstPtr(); // interior pointer
    Lucee::FieldPtr<T> gPtr = this->createPtr(); // ghost pointer

    int idx[NDIM];
// loop over region and copy to ghost layer of cells
    Lucee::RowMajorSequencer<NDIM> loSeq(lowerRgn);
    while (loSeq.step())
    {
      loSeq.fillWithIndex(idx);
      this->setPtr(ptr, idx); // location to copy from
// bump index in 'dir' to get to correction location to copy into
      idx[dir] += stride;
      this->setPtr(gPtr, idx);

// copy data
      for (unsigned k=0; k<numComponents; ++k)
        gPtr[k] = ptr[k];
    }

// create a region to represent upper interior layer of cells
    for (unsigned i=0; i<NDIM; ++i)
    { // whole interior region
      lo[i] = getLower(i);
      up[i] = getUpper(i);
    }
// adjust cells along 'dir' direction so it represents region to be copied
    lo[dir] = up[dir]-lowerGhost[dir];
    Lucee::Region<NDIM, int> upperRgn(lo, up);

// loop over region and copy to ghost layer of cells
    Lucee::RowMajorSequencer<NDIM> upSeq(upperRgn);
    while (upSeq.step())
    {
      upSeq.fillWithIndex(idx);
      this->setPtr(ptr, idx); // location to copy from
// bump index in 'dir' to get to correction location to copy into
      idx[dir] -= stride;
      this->setPtr(gPtr, idx);

// copy data
      for (unsigned k=0; k<numComponents; ++k)
        gPtr[k] = ptr[k];
    }

    return *this;
  }

  template <unsigned NDIM, typename T>
  Field<NDIM, T>&
  Field<NDIM, T>::applyCopyBc(unsigned dir, unsigned side)
  {
    int lo[NDIM], up[NDIM];

// create a region to represent the ghost layer
    for (unsigned i=0; i<NDIM; ++i)
    { // whole region, including extended region
      lo[i] = getGlobalLowerExt(i);
      up[i] = getGlobalUpperExt(i);
    }
// adjust region so it only indexes the ghost cells
    if (side == Lucee::LOWER_SIDE)
    { // lower side
      up[dir] = getGlobalLower(dir);
    }
    else
    { // upper side
      lo[dir] = getGlobalUpper(dir);
    }
// region must be local to processor
    Lucee::Region<NDIM, int> gstRgn = getExtRegion().intersect(
      Lucee::Region<NDIM, int>(lo, up));

    Lucee::ConstFieldPtr<T> ptr = this->createConstPtr(); // interior pointer
    Lucee::FieldPtr<T> gPtr = this->createPtr(); // ghost pointer

    int idx[NDIM];
// loop over region and from first layer of cells
    Lucee::RowMajorSequencer<NDIM> seq(gstRgn);
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      this->setPtr(gPtr, idx); // index into ghost cell
      if (side == 0)
      { // lower
        idx[dir] = getLower(dir); // first location in domain
      }
      else
      { // upper
        idx[dir] = getUpper(dir)-1; // last location in domain
      }
      this->setPtr(ptr, idx);
      for (unsigned k=0; k<numComponents; ++k)
        gPtr[k] = ptr[k];
    }

    return *this;
  }

  template <unsigned NDIM, typename T>
  void
  Field<NDIM, T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("clear", luaClear);
    lfm.appendFunc("scale", luaScale);
    lfm.appendFunc("copy", luaCopy);
    lfm.appendFunc("accumulate", luaAccumulate);
    lfm.appendFunc("combine", luaCombine);
    lfm.appendFunc("hasNan", luaHasNan);
    lfm.appendFunc("applyPeriodicBc", luaApplyPeriodicBc);
    lfm.appendFunc("applyCopyBc", luaApplyCopyBc);
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
  Field<NDIM, T>::luaScale(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (! lua_isnumber(L, 2))
    {
      Lucee::Except lce("Field::luaScale: Must provide a number to 'scale' method");
      throw lce;
    }
    T num = (T) lua_tonumber(L, 2);
    (*fld) *= num;
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
// number of arguments in call
    unsigned nArgs = lua_gettop(L);

    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);

// accumulate each field into this one
    for (int nF = 0; nF < (nArgs-1)/2; ++nF)
    {
      if (! lua_isnumber(L, 2*nF+2))
      {
        Lucee::Except lce("Field::luaAccumulate: Must provide a number to 'accumulate' method");
        throw lce;
      }
      T coeff = (T) lua_tonumber(L, 2*nF+2); // coeff for accumulation

      if (lua_type(L, 2*nF+3) != LUA_TUSERDATA)
      {
        Lucee::Except lce("Field::luaAccumulate: Must provide a field to 'accumulate' method");
        throw lce;
      }
      Lucee::PointerHolder<Field<NDIM, T> > *fldPtr =
        (Lucee::PointerHolder<Field<NDIM, T> >*) lua_touserdata(L, 2*nF+3); // field to accumulate
      fld->accumulate(coeff, *fldPtr->pointer);
    }

    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaCombine(lua_State *L)
  {
// number of arguments in call
    unsigned nArgs = lua_gettop(L);

    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);

// clear out field
    (*fld) = (T) 0;

// accumulate each field into this one
    for (int nF = 0; nF < (nArgs-1)/2; ++nF)
    {
      if (! lua_isnumber(L, 2*nF+2))
      {
        Lucee::Except lce("Field::luaAccumulate: Must provide a number to 'combine' method");
        throw lce;
      }
      T coeff = (T) lua_tonumber(L, 2*nF+2); // coeff for accumulation

      if (lua_type(L, 2*nF+3) != LUA_TUSERDATA)
      {
        Lucee::Except lce("Field::luaAccumulate: Must provide a field to 'combine' method");
        throw lce;
      }
      Lucee::PointerHolder<Field<NDIM, T> > *fldPtr =
        (Lucee::PointerHolder<Field<NDIM, T> >*) lua_touserdata(L, 2*nF+3); // field to accumulate
      fld->accumulate(coeff, *fldPtr->pointer);
    }

    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaHasNan(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);

    int hasNan = false;
// loop over field, checking for nans
    Lucee::ConstFieldPtr<T> ptr = fld->createConstPtr();
    Lucee::RowMajorSequencer<NDIM> seq(fld->getExtRegion());
    while (seq.step())
    {
      fld->setPtr(ptr, seq.getIndex());
      for (unsigned k=0; k<ptr.getNumComponents(); ++k)
      {
        if (Lucee::isNan(ptr[k])) 
        {
          hasNan = true;
          break;
        }
      }
    }

// make sure all processors return same answer
    TxCommBase *comm = fld->getComm();
    int globalResult;
    comm->allreduce(1, &hasNan, &globalResult, TX_AND);

// push results on stack
    if (globalResult)
      lua_pushboolean(L, 1);
    else
      lua_pushboolean(L, 0);

    return 1;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaApplyPeriodicBc(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (! lua_isnumber(L, 2))
    {
      Lucee::Except lce("Field::luaApplyPeriodicBc: Must provide a number to 'applyPeriodicBc' method");
      throw lce;
    }
// determine direction in which to apply periodic BCs
    int dir = (int) lua_tonumber(L, 2);
    if (dir<0 || dir >= NDIM)
    { // incorrect direction specified
      Lucee::Except lce("Field::luaApplyPeriodicBc: Direction must be one of ");
      for (unsigned i=0; i<NDIM-1; ++i)
        lce << i << ", ";
      lce << NDIM-1 << ".";
      lce << " '" << dir << "' specified instead" << std::endl;
      throw lce;      
    }
// apply boundary conditions
    fld->applyPeriodicBc(dir);

    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  Field<NDIM, T>::luaApplyCopyBc(lua_State *L)
  {
    Field<NDIM, T> *fld
      = Lucee::PointerHolder<Field<NDIM, T> >::getObj(L);
    if (! lua_isnumber(L, 2))
    {
      Lucee::Except lce("Field::luaApplyCopyBc: Must provide a number to 'applyCopyBc' method");
      throw lce;
    }
// determine direction in which to apply copy BCs
    int dir = (int) lua_tonumber(L, 2);
    if (dir<0 || dir >= NDIM)
    { // incorrect direction specified
      Lucee::Except lce("Field::luaApplyCopyBc: Direction must be one of ");
      for (unsigned i=0; i<NDIM-1; ++i)
        lce << i << ", ";
      lce << NDIM-1 << ".";
      lce << " '" << dir << "' specified instead" << std::endl;
      throw lce;      
    }

    if (! lua_isstring(L, 3))
    {
      Lucee::Except lce("Field::luaApplyCopyBc: Must provide a side to 'applyCopyBc' method.");
      lce << " Should be one of 'lower' or 'upper'." << std::endl;
      throw lce;
    }
// determine side in which to apply periodic BCs
    std::string ss = lua_tostring(L, 3);
    unsigned side = 1;
    if (ss == "lower")
      side = Lucee::LOWER_SIDE;
    else if (ss == "upper")
      side = Lucee::UPPER_SIDE;
    else
      throw Lucee::Except("Field::luaApplyCopyBc: side should be one of \"lower\" or \"upper\".");
// apply boundary conditions
    fld->applyCopyBc(dir, side);

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
