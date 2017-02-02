/**
 * @file	LcStructGridField.cpp
 *
 * @brief	StructGridFields are fields that live on structured grids.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathLib.h>
#include <LcPointerHolder.h>
#include <LcStructGridField.h>
#include <LcMatrix.h>

// loki includes
#include <loki/Singleton.h>

// boost includes
#include <boost/shared_array.hpp>

// std includes
#include <fstream>

namespace Lucee
{
// flags for data location
  static const unsigned CELL_CENTER_LOC = 0;
  static const unsigned VERTEX_LOC = 1;

// names used in registration system
  template <> const char *StructGridField<1, double>::id = "Field1D";
  template <> const char *StructGridField<2, double>::id = "Field2D";
  template <> const char *StructGridField<3, double>::id = "Field3D";
  template <> const char *StructGridField<4, double>::id = "Field4D";
  template <> const char *StructGridField<5, double>::id = "Field5D";

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField()
    : Lucee::Field<NDIM, T>(), grid(0)
  {
    this->setName(StructGridField<NDIM, double>::id);
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(Lucee::StructuredGridBase<NDIM>* grid, unsigned nc,
        int lg[NDIM], int ug[NDIM])
    : Lucee::Field<NDIM, T>(grid->getGlobalRegion(), grid->getLocalRegion(), nc, lg, ug, (T) 0),
      grid(grid)
  {
    this->setName(StructGridField<NDIM, double>::id);
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(const StructGridField<NDIM, T>& fld)
    : Lucee::Field<NDIM, T>(fld), grid(fld.grid)
  {
    this->setName(StructGridField<NDIM, double>::id);
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>&
  StructGridField<NDIM, T>::operator=(const StructGridField<NDIM, T>& fld)
  {
    if (&fld == this)
      return *this;
// call base class assignment operator
    Field<NDIM, T>::operator=(fld);
    grid = fld.grid;

    return *this;
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>&
  StructGridField<NDIM, T>::operator=(const T& val)
  {
    Field<NDIM, T>::operator=(val);

    return *this;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::divergence(Lucee::StructGridField<NDIM, T>& div) const
  {
// WARNING: THIS CODE PRESENTLY ONLY WORKS FOR RECTCART GRIDS

    if (this->getNumComponents() != NDIM)
    {
      Lucee::Except lce("StructGridField::divergence: Incorrect number of components. Should be ");
      lce << NDIM << " but has " << this->getNumComponents() << " components." ;
      throw lce;
    }

// create various iterators
    Lucee::FieldPtr<T> divPtr = div.createPtr();
    Lucee::ConstFieldPtr<T> rPtr = this->createConstPtr();
    Lucee::ConstFieldPtr<T> lPtr = this->createConstPtr();

    int idx[NDIM]; // for indexing

    div = 0.0; // clear divergence field
    for (unsigned n=0; n<NDIM; ++n)
    {
      double dx1 = 0.5/grid->getDx(n);
      Lucee::RowMajorSequencer<NDIM> seq(this->getRegion());
      while (seq.step())
      {
        seq.fillWithIndex(idx);
        div.setPtr(divPtr, idx); // current cell
        idx[n] = idx[n]+1;
        this->setPtr(rPtr, idx); // right cell
        idx[n] = idx[n]-1-1;
        this->setPtr(lPtr, idx); // left cell
        divPtr[0] += dx1*(rPtr[n]-lPtr[n]);
      }
    }

  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::copyFromCptr(T *cptr)
  {
    Lucee::FieldPtr<T> ptr = this->createPtr(); // pointer to help in setting field
    Lucee::Region<NDIM, int> extRgn = this->getExtRegion(); // loop over extended region
    Lucee::RowMajorSequencer<NDIM> seq(extRgn);
    int idx[NDIM];
    double xc[3];
    unsigned numOut = this->getNumComponents();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      this->setPtr(ptr, idx);
      grid->setIndex(idx);
      if (dataLoc == VERTEX_LOC)
        grid->getVertex(xc); // vertex coordinate
      else
        grid->getCentroid(xc); // cell center coordinate
      for (int i=0; i<numOut; ++i)
      {
        ptr[i] = cptr[i];
      }
      cptr += numOut;
    }
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::copyToCptr(T *cptr)
  {
    Lucee::FieldPtr<T> ptr = this->createPtr(); // pointer to help in setting field
    Lucee::Region<NDIM, int> extRgn = this->getExtRegion(); // loop over extended region
    Lucee::RowMajorSequencer<NDIM> seq(extRgn);
    int idx[NDIM];
    double xc[3];
    unsigned numOut = this->getNumComponents();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      this->setPtr(ptr, idx);
      grid->setIndex(idx);
      if (dataLoc == VERTEX_LOC)
        grid->getVertex(xc); // vertex coordinate
      else
        grid->getCentroid(xc); // cell center coordinate
      for (int i=0; i<numOut; ++i)
      {
        cptr[i] = ptr[i];
      }
      cptr += numOut;
    }
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
// get name of grid on which field lives
    if (tbl.template hasObject<Lucee::StructuredGridBase<NDIM> >("onGrid"))
      grid = &tbl.template
        getObject<Lucee::StructuredGridBase<NDIM> > ("onGrid");
    else
      throw Lucee::Except("StructGridField::readInput: must specify 'onGrid', the grid on which field lives");

// set communicators and I/O flag
    this->setComm(grid->getComm()); this->setMomComm(grid->getMomComm());
    this->setIsSafeToWrite(grid->isSafeToWrite());

// check where data should be located
    dataLoc = CELL_CENTER_LOC;
    if (tbl.hasString("location"))
    {
      if (tbl.getString("location") == "vertex")
        dataLoc = VERTEX_LOC;
      else if (tbl.getString("location") == "center")
        dataLoc = CELL_CENTER_LOC;
      else
      {
        Lucee::Except lce(
          "StructGridField::readInput: 'location' must be one of 'vertex' or 'center'. Provided ");
        lce << tbl.getString("location") << " instead" << std::endl;
        throw lce;
      }
    }

    unsigned numComponents = 1;
// get number of components in field
    if (tbl.hasNumber("numComponents"))
      numComponents = tbl.getNumber("numComponents");

    int lowerGhost[NDIM], upperGhost[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
// read number of ghost cells
    if (tbl.hasNumVec("ghost"))
    {
      std::vector<double> gstDbl = tbl.getNumVec("ghost");
      if (gstDbl.size() != 2)
        throw Lucee::Except("StructGridField::readInput: must specify exactly two entries in 'ghost' table");
      for (unsigned i=0; i<NDIM; ++i)
      {
        lowerGhost[i] = gstDbl[0];
        upperGhost[i] = gstDbl[1];
      }
    }

    for (unsigned i=0; i<NDIM; ++i)
    {
      lowerWriteGhost[i] = 0;
      upperWriteGhost[i] = 0;
    }
// read number of ghost cells to write
    if (tbl.hasNumVec("writeGhost"))
    {
      std::vector<double> gstDbl = tbl.getNumVec("writeGhost");
      if (gstDbl.size() != 2)
        throw Lucee::Except("StructGridField::readInput: must specify exactly two entries in 'writeGhost' table");
      for (unsigned i=0; i<NDIM; ++i)
      {
        lowerWriteGhost[i] = gstDbl[0];
        upperWriteGhost[i] = gstDbl[1];
      }
    }

// check if this field lives in parallel
    bool isPar = true; // by default alway parallel
    if (tbl.hasBool("decompose"))
      isPar = tbl.getBool("decompose");

// global region for field
    typename Lucee::Region<NDIM, int> globalRgn = grid->getGlobalRegion();
// local region for field
    typename Lucee::Region<NDIM, int> localRgn = globalRgn;
    if (isPar)
      localRgn = grid->getLocalRegion();
// create new field and copy data to self
    Field<NDIM, T>::operator=(
      Field<NDIM, T>(globalRgn, localRgn, numComponents, lowerGhost, upperGhost));
  }

  template <unsigned NDIM, typename T>
  TxIoNodeType
  StructGridField<NDIM, T>::writeToFile(TxIoBase& io, TxIoNodeType& node, const std::string& nm)
  {
// first write the field data to file
    TxIoNodeType dn 
      = Lucee::Field<NDIM, T>::writeToFileWithGhost(io, node, 
        lowerWriteGhost, upperWriteGhost, "StructGridField");
// annotate with viz-schema marks
    io.writeAttribute(dn, "vsType", "variable");
    io.writeAttribute(dn, "vsMesh", "StructGrid");
    if (dataLoc == CELL_CENTER_LOC)
      io.writeAttribute(dn, "vsCentering", "zonal");
    else
      io.writeAttribute(dn, "vsCentering", "nodal");
// now write out grid
    grid->writeToFile(io, node, "StructGrid");

    return node;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::writeToTxtFile(std::ofstream& txtFl)
  {
// create sequencer to loop over complete region
    Lucee::RowMajorSequencer<NDIM> seq(this->getRegion());
    
    double xc[3]; // for cell-center coordinates
    int idx[NDIM]; // for indexing
    Lucee::ConstFieldPtr<T> rPtr = this->createConstPtr(); // pointer to data location
    while (seq.step())
    {
      seq.fillWithIndex(idx);
// set pointers in data array and grid
      this->setPtr(rPtr, idx);
      grid->setIndex(idx);
// get cell-center coordinates
      grid->getCentroid(xc);

// write coordinates of cell-center
      for (unsigned i=0; i<NDIM; ++i)
        txtFl << xc[i] << " ";
// now write out actual data at this location
      for (unsigned i=0; i<this->getNumComponents(); ++i)
        txtFl << rPtr[i] << " ";
      txtFl << std::endl;
    }
  }

  template <unsigned NDIM, typename T>
  TxIoNodeType
  StructGridField<NDIM, T>::readFromFile(TxIoBase& io, TxIoNodeType& node,
    const std::string& nm)
  {
    TxIoNodeType dn 
      = Lucee::Field<NDIM, T>::readFromFile(io, node, "StructGridField");
    return dn;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::sync()
  {
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    int rank = comm->getRank(); // current rank
    int lg[NDIM], ug[NDIM];
    this->fillWithGhosts(lg, ug); // ghost cells

// initiate (non-blocking) receive of data from our neighbors
    startRecv(rank, lg, ug);
// send data over
    send(rank, lg, ug);
// finish receiving data as recieve was non-blocking
    finishRecv();
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::startRecv(unsigned rank, int lg[NDIM], int ug[NDIM])
  {
// finish receiving if we are already receiving
    if (isReceiving) finishRecv();

// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();

    Lucee::Region<NDIM, int> localBox = this->getRegion();
// get neighbor information
    std::vector<unsigned> neigh = grid->getRecvNeighbors(rank, lg, ug);
    std::vector<unsigned>::const_iterator niter;
// receive stuff from our neighbors
    for (niter = neigh.begin(); niter != neigh.end(); ++niter) 
    {
      int srcrank = grid->getRgnRank(*niter);
      Lucee::Region<NDIM, int> otherBox = grid->getNeighborRgn(*niter);
      int tag = 1000000;
      if (rank == srcrank) 
//  skip if we're receiving from ourselves
        continue;

      if (srcrank != *niter) 
      { // periodic direction
//  receiving a wrapped boundary block, tag according to direction.
//  receiving from the negative side = -1, receiving from the positive
//  side = +1
        int tagmod = 1;
        for (unsigned i=0; i<NDIM; ++i)
        {
          if (otherBox.getLower(i) > localBox.getLower(i)) 
            tag += tagmod;
          else if (otherBox.getLower(i) < localBox.getLower(i)) 
            tag -= tagmod;
          tagmod *= 10;
        }
      }
//  set up a communication
      Lucee::Region<NDIM, int> recvRgn = otherBox.intersect( this->getExtRegion() );
      if ( recvRgn.getVolume() == 0 )
        continue; // nothing do if no intersection

      TxMsgStatus ms = comm->template
        startRecv<T>(this->getNumComponents()*recvRgn.getVolume(), srcrank, tag);
      msgStatus[*niter] = ms;  // key in msgStatus map is region we expect receive
    }
    isReceiving = true;
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::send(unsigned rank, int lg[NDIM], int ug[NDIM])
  {
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();

    Lucee::Region<NDIM, int> localBox = this->getRegion();

    std::vector<unsigned> neigh = grid->getSendNeighbors(rank, lg, ug);
    std::vector<unsigned>::const_iterator niter;
    for (niter = neigh.begin(); niter != neigh.end(); ++niter)  
    {
      int destrank = grid->getRgnRank(*niter);
      Lucee::Region<NDIM, int> otherBox = grid->getNeighborRgn(*niter);
      int tag = 1000000;

      if (rank == destrank) 
      {
//  destination is the same as the source.  Just copy it over
        Lucee::Region<NDIM, int> fakeBox = grid->getNeighborRgn(*niter);
        int lCor[NDIM], uCor[NDIM];
        for (unsigned i=0; i<NDIM; ++i) 
        {
          lCor[i] = fakeBox.getLower(i) - localBox.getLower(i);
          uCor[i] = localBox.getLower(i) - fakeBox.getLower(i);
        }
        Lucee::Region<NDIM, int> sendRgn =  localBox.intersect(fakeBox.extend(lg, ug));
        Lucee::Region<NDIM, int> recvRgn = (localBox.intersect(
            fakeBox.extend(lg, ug))).extend(lCor, uCor);

        Lucee::RowMajorSequencer<NDIM> sendSeq(sendRgn);
        Lucee::RowMajorSequencer<NDIM> recvSeq(recvRgn);
        Lucee::FieldPtr<T> sendItr = this->createPtr();
        Lucee::FieldPtr<T> recvItr = this->createPtr();
        while ( sendSeq.step() && recvSeq.step() ) 
        {
          this->setPtr(sendItr, sendSeq.getIndex());
          this->setPtr(recvItr, recvSeq.getIndex());
          for (unsigned k=0; k<this->getNumComponents(); ++k)
            recvItr[k] = sendItr[k];
        }
        continue;
      }
      if (destrank != *niter) 
      {
// sending a wrapped boundary block, tag according to direction.
// sending to the negative side = +1, sending to the positive side =
// -1
        int tagmod = 1;
        for (unsigned i=0; i<NDIM; ++i) 
        {
          if (otherBox.getLower(i) > localBox.getLower(i))
            tag -= tagmod;
          else if (otherBox.getLower(i) < localBox.getLower(i))
            tag += tagmod;
          tagmod *= 10;
        }
      }
      Lucee::Region<NDIM, int> sendRgn = localBox.intersect(
        grid->getNeighborRgn(*niter).extend(lg, ug));
      if ( sendRgn.getVolume() == 0 )
        continue;
      std::vector<T> sendVec(this->getNumComponents()*sendRgn.getVolume());
      Lucee::RowMajorSequencer<NDIM> seq(sendRgn);
      int i = 0;
      Lucee::FieldPtr<T> itr = this->createPtr();
      while ( seq.step() ) 
      {
        this->setPtr(itr, seq.getIndex());
        for (unsigned k=0; k<this->getNumComponents(); ++k)
          sendVec[i++] = itr[k];
      }
      comm->template
        send<T>(this->getNumComponents()*sendRgn.getVolume(), &sendVec[0], destrank, tag);
    }
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::finishRecv()
  {
    if (!isReceiving) return; // nothing to receive

// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();

// complete receives and copy data into ghost cells
    std::map<int, TxMsgStatus>::iterator i = msgStatus.begin();
    for ( ; i != msgStatus.end(); ++i)
    {
      Lucee::Region<NDIM, int> recvRgn = this->getExtRegion().intersect(
        grid->getNeighborRgn(i->first)); // ghost cell region
      T* recvData = (T*) comm->finishRecv(i->second);

// copy data over: this works because send() function uses same
// ordering of data as we use below.
      Lucee::RowMajorSequencer<NDIM> seq(recvRgn);
      Lucee::FieldPtr<T> ptr = this->createPtr();
      unsigned loc = 0;
      while ( seq.step() )
      {
        this->setPtr(ptr, seq.getIndex());
        for (unsigned k=0; k<this->getNumComponents(); ++k)
          ptr[k] = recvData[loc++];
      }
      delete [] recvData;
    }
    msgStatus.erase(msgStatus.begin(), msgStatus.end());
    isReceiving = false; // we are done
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    Field<NDIM, T>::appendLuaCallableMethods(lfm);
// now append local methods
    lfm.appendFunc("set", luaSet);
    lfm.appendFunc("alias", luaAlias);
    lfm.appendFunc("duplicate", luaDuplicate);
    lfm.appendFunc("div", luaDivergence);
    lfm.appendFunc("applyFuncBc", luaSetGhost);
    lfm.appendFunc("copy_from_cptr", luaCopyFromCptr);
    lfm.appendFunc("copy_to_cptr", luaCopyToCptr);
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaSet(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);

    if (! lua_isfunction(L, 2))
    {
      Lucee::Except lce("StructGridField::luaSet: Must provide a Lua function to 'set' method");
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_pop(L, 1); // WHY IS THIS NEEDED? THIS LOOKS LIKE A POTENTIAL STACK CORRUPTION

    sgf->setFromLuaFunction(L, fnRef);
    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaSetGhost(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);

    if (! lua_isnumber(L, 3))
    {
      Lucee::Except lce(
        "StructGridField::luaSetGhost: Must provide a number to 'applyFuncBc' method");
      throw lce;
    }

// determine direction in which to apply function BCs
    int dir = (int) lua_tonumber(L, 3);
    if (dir<0 || dir >= NDIM)
    { // incorrect direction specified
      Lucee::Except lce(
        "StructGridField::luaSetGhost: Direction must be one of ");
      for (unsigned i=0; i<NDIM-1; ++i)
        lce << i << ", ";
      lce << NDIM-1 << ".";
      lce << " '" << dir << "' specified instead" << std::endl;
      throw lce;      
    }

    if (! lua_isstring(L, 4))
    {
      Lucee::Except lce("StructGridField::luaFuncBc: Must provide a side to 'applyFuncBc' method.");
      lce << " Should be one of 'lower' or 'upper'." << std::endl;
      throw lce;
    }
// determine side in which to apply periodic BCs
    std::string ss = lua_tostring(L, 3);
    unsigned side = 1;
    if (ss == "lower")
      side = 0;
    else if (ss == "upper")
      side = 1;
    else
      throw Lucee::Except("StructGridField::luaFuncBc: side should be one of \"lower\" or \"upper\".");

    if (! lua_isfunction(L, 2))
    {
      Lucee::Except lce(
        "StructGridField::luaSetGhost: Must provide a Lua function to 'applyFuncBc' method");
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);

    sgf->setGhostFromLuaFunction(L, fnRef, dir, side);
    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaAlias(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);
    if (! lua_isnumber(L, 2) && ! lua_isnumber(L, 3))
    {
      Lucee::Except lce("StructGridField::luaAlias: 'alias' method must be passed half-open interval as [s,e)");
      throw lce;
    }
    unsigned s = lua_tonumber(L, 2);
    unsigned e = lua_tonumber(L, 3);

    Lucee::Field<NDIM, T> aliasFld = sgf->getSubCompView(s, e);
    Lucee::StructGridField<NDIM, T>* aliasSgf = 
      new Lucee::StructGridField<NDIM, T>(aliasFld, *sgf->grid);
// copy data for ghost I/O
    for (unsigned d=0; d<NDIM; ++d)
    {
      aliasSgf->lowerWriteGhost[d] = sgf->lowerWriteGhost[d];
      aliasSgf->upperWriteGhost[d] = sgf->upperWriteGhost[d];
    }

    size_t nbytes = sizeof(Lucee::PointerHolder<StructGridField<NDIM, T> >);
    Lucee::PointerHolder<StructGridField<NDIM, T> > *ph =
      (Lucee::PointerHolder<StructGridField<NDIM, T> >*) lua_newuserdata(L, nbytes);
    ph->pointer = aliasSgf;
    ph->pointer->template setBaseType<Lucee::DataStructIfc>();
    ph->pointer->template setDerivedType<StructGridField<NDIM, T> >();

    luaL_getmetatable(L, typeid(StructGridField<NDIM, T>).name());
    lua_setmetatable(L, -2);

    return 1;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaDuplicate(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);

    Lucee::StructGridField<NDIM, T>* aliasSgf = 
      new Lucee::StructGridField<NDIM, T>(
        sgf->duplicate(), *sgf->grid);

    aliasSgf->dataLoc = sgf->dataLoc;
    for (unsigned d=0; d<NDIM; ++d)
    {
      aliasSgf->lowerWriteGhost[d] = sgf->lowerWriteGhost[d];
      aliasSgf->upperWriteGhost[d] = sgf->upperWriteGhost[d];
    }

    size_t nbytes = sizeof(Lucee::PointerHolder<StructGridField<NDIM, T> >);
    Lucee::PointerHolder<StructGridField<NDIM, T> > *ph =
      (Lucee::PointerHolder<StructGridField<NDIM, T> >*) lua_newuserdata(L, nbytes);
    ph->pointer = aliasSgf;
    ph->pointer->template setBaseType<Lucee::DataStructIfc>();
    ph->pointer->template setDerivedType<StructGridField<NDIM, T> >();

    luaL_getmetatable(L, typeid(StructGridField<NDIM, T>).name());
    lua_setmetatable(L, -2);

    return 1;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaDivergence(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);
    if (lua_type(L, 2) != LUA_TUSERDATA)
    {
      Lucee::Except lce("StructGridField::luaDivergence: Must provide a field to 'div' method");
      throw lce;
    }
    Lucee::PointerHolder<StructGridField<NDIM, T> > *fldPtr =
      (Lucee::PointerHolder<StructGridField<NDIM, T> >*) lua_touserdata(L, 2);
// compute divergence
    sgf->divergence(*fldPtr->pointer);

    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaCopyFromCptr(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);

    if (!lua_isuserdata(L, 2)) {
      Lucee::Except lce(
        "StructGridField::luaCopyFromCptr: Must provide userdata");
      throw lce;
    }
    
    sgf->copyFromCptr((T*) lua_touserdata(L, 2));
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaCopyToCptr(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);

    if (!lua_isuserdata(L, 2)) {
      Lucee::Except lce(
        "StructGridField::luaCopyFromCptr: Must provide userdata");
      throw lce;
    }
    
    sgf->copyToCptr((T*) lua_touserdata(L, 2));
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::setFromLuaFunction(lua_State *L, int ref)
  {
    Lucee::FieldPtr<T> ptr = this->createPtr(); // pointer to help in setting field
    Lucee::Region<NDIM, int> extRgn = this->getExtRegion(); // loop over extended region
    Lucee::RowMajorSequencer<NDIM> seq(extRgn);
    int idx[NDIM];
    double xc[3];
    unsigned numOut = this->getNumComponents();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      this->setPtr(ptr, idx);
      grid->setIndex(idx);
      if (dataLoc == VERTEX_LOC)
        grid->getVertex(xc); // vertex coordinate
      else
        grid->getCentroid(xc); // cell center coordinate
// push function object on stack
      lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// ensure enough space on stack
      lua_checkstack(L, numOut+3);
// push variables on stack
      for (unsigned i=0; i<3; ++i){
        lua_pushnumber(L, xc[i]);
      }
      for (int i=0; i<numOut; ++i)
      {
        lua_pushnumber(L, ptr[i]);
      }
      if (lua_pcall(L, 3+numOut, numOut, 0) != 0)
      {
        Lucee::Except lce("StructGridField::setFromLuaFunction: ");
        lce << "Problem evaluating function supplied to 'set' method"
            << std::endl;
        std::string err(lua_tostring(L, -1));
        lua_pop(L, 1);
        lce << "[" << err << "]";
        throw lce;
      }
// fetch results
      for (int i=-numOut; i<0; ++i)
      {
        if (!lua_isnumber(L, i))
          throw Lucee::Except("StructGridField::setFromLuaFunction: Return value not a number");
        ptr[numOut+i] = lua_tonumber(L, i);
      }
      lua_pop(L, 1);
    }
  }

  template <unsigned NDIM, typename T>
  void
  StructGridField<NDIM, T>::setGhostFromLuaFunction(lua_State *L, int ref,
    unsigned dir, unsigned side)
  {
    Lucee::FieldPtr<T> ptr = this->createPtr(); // pointer to help in setting field
    int lo[NDIM], up[NDIM];

// create a region to represent ghost layer
    for (unsigned i=0; i<NDIM; ++i)
    { // whole region, including extended region
      lo[i] = this->getGlobalLowerExt(i);
      up[i] = this->getGlobalUpperExt(i);
    }
// adjust region so it only indexes ghost cells
    if (side == 0)
    { // lower side
      up[dir] = this->getGlobalLower(dir);
    }
    else
    { // upper side
      lo[dir] = this->getGlobalUpper(dir);
    }
// region must be local to processor
    Lucee::Region<NDIM, int> gstRgn = this->getExtRegion().intersect(
      Lucee::Region<NDIM, int>(lo, up));

    Lucee::RowMajorSequencer<NDIM> seq(gstRgn);
    int idx[NDIM];
    double xc[3];
    unsigned numOut = this->getNumComponents();
    while (seq.step())
    {
      seq.fillWithIndex(idx);
      this->setPtr(ptr, idx);
      grid->setIndex(idx);
      grid->getCentroid(xc); // cell center coordinate
// push function object on stack
      lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push variables on stack
      for (unsigned i=0; i<3; ++i)
        lua_pushnumber(L, xc[i]);
      if (lua_pcall(L, 3, numOut, 0) != 0)
      {
        Lucee::Except lce("StructGridField::setGhostFromLuaFunction: ");
        lce << "Problem evaluating function supplied to 'set' method";
        throw lce;
      }
// fetch results
      for (int i=-numOut; i<0; ++i)
      {
        if (!lua_isnumber(L, i))
          throw Lucee::Except("StructGridField::setGhostFromLuaFunction: Return value not a number");
        ptr[numOut+i] = lua_tonumber(L, i);
      }
      lua_pop(L, 1);
    }
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(const Field<NDIM, T>& fld, Lucee::StructuredGridBase<NDIM>& grd)
    : Lucee::Field<NDIM, T>(fld), grid(&grd)
  {
    this->setName(StructGridField<NDIM, double>::id);
  }

// instantiations
  template class StructGridField<1, int>;
  template class StructGridField<2, int>;
  template class StructGridField<3, int>;
  template class StructGridField<4, int>;
  template class StructGridField<5, int>;

  template class StructGridField<1, float>;
  template class StructGridField<2, float>;
  template class StructGridField<3, float>;
  template class StructGridField<4, float>;
  template class StructGridField<5, float>;

  template class StructGridField<1, double>;
  template class StructGridField<2, double>;
  template class StructGridField<3, double>;
  template class StructGridField<4, double>;
  template class StructGridField<5, double>;

}
