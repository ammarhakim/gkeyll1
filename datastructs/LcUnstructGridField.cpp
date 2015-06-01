/**
 * @file	LcUnstructGridField.cpp
 *
 * @brief	UnstructGridFields are fields that live on structured grids.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcMathLib.h>
#include <LcPointerHolder.h>
#include <LcUnstructGridField.h>
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
  template<> const char *UnstructGridField<1, double>::id = "UField1D";
  template<> const char *UnstructGridField<2, double>::id = "UField2D";
  template<> const char *UnstructGridField<3, double>::id = "UField3D";

  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>::UnstructGridField() :
      Lucee::Field<NDIM, T>(), grid(0)
  {
    this->setName(UnstructGridField<NDIM, double>::id);
  }

  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>::UnstructGridField(
      Lucee::UnstructuredGrid<NDIM>* grid, unsigned nc) : grid(grid)
  {
    this->setName(UnstructGridField<NDIM, double>::id);
  }

  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>::UnstructGridField(
      const UnstructGridField<NDIM, T>& fld) :
      Lucee::Field<NDIM, T>(fld), grid(fld.grid)
  {
    this->setName(UnstructGridField<NDIM, double>::id);
  }

  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>&
  UnstructGridField<NDIM, T>::operator=(const UnstructGridField<NDIM, T>& fld)
  {
    if (&fld == this)
      return *this;
// call base class assignment operator
    Field<NDIM, T>::operator=(fld);
    grid = fld.grid;

    return *this;
  }

  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>&
  UnstructGridField<NDIM, T>::operator=(const T& val)
  {
    Field<NDIM, T>::operator=(val);

    return *this;
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
// get name of grid on which field lives
    if (tbl.template hasObject<Lucee::UnstructuredGrid<NDIM> >("onGrid"))
      grid = &tbl.template getObject<Lucee::UnstructuredGrid<NDIM> >("onGrid");
    else
      throw Lucee::Except(
          "UnstructGridField::readInput: must specify 'onGrid', the grid on which field lives");

// set communicators and I/O flag
    this->setComm(grid->getComm());
    this->setMomComm(grid->getMomComm());
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
            "UnstructGridField::readInput: 'location' must be one of 'vertex' or 'center'. Provided ");
        lce << tbl.getString("location") << " instead" << std::endl;
        throw lce;
      }
    }

    unsigned numComponents = 1;
// get number of components in field
    if (tbl.hasNumber("numComponents"))
      numComponents = tbl.getNumber("numComponents");

    int lowerGhost[NDIM], upperGhost[NDIM];
    for (unsigned i = 0; i < NDIM; ++i)
    {
      lowerGhost[i] = 0;
      upperGhost[i] = 0;
    }
// read number of ghost cells
    if (tbl.hasNumVec("ghost"))
    {
      std::vector<double> gstDbl = tbl.getNumVec("ghost");
      if (gstDbl.size() != 2)
        throw Lucee::Except(
            "UnstructGridField::readInput: must specify exactly two entries in 'ghost' table");
      for (unsigned i = 0; i < NDIM; ++i)
      {
        lowerGhost[i] = gstDbl[0];
        upperGhost[i] = gstDbl[1];
      }
    }

// read number of ghost cells to write
    if (tbl.hasNumVec("writeGhost"))
    {
      std::vector<double> gstDbl = tbl.getNumVec("writeGhost");
      if (gstDbl.size() != 2)
        throw Lucee::Except(
            "UnstructGridField::readInput: must specify exactly two entries in 'writeGhost' table");

    }

// check if this field lives in parallel
    bool isPar = true; // by default alway parallel
    if (tbl.hasBool("decompose"))
      isPar = tbl.getBool("decompose");

    //TODO: Fields are never set up

#if 0

// global region for field
    typename Lucee::Region<NDIM, int> globalRgn = grid->getGlobalRegion();
// local region for field
    typename Lucee::Region<NDIM, int> localRgn = globalRgn;
    if (isPar)
      localRgn = grid->getLocalRegion();
// create new field and copy data to self
    Field<NDIM, T>::operator=(
        Field<NDIM, T>(globalRgn, localRgn, numComponents));
#endif
  }

  template<unsigned NDIM, typename T>
  TxIoNodeType UnstructGridField<NDIM, T>::writeToFile(TxIoBase& io,
      TxIoNodeType& node, const std::string& nm)
  {
// first write the field data to file
    TxIoNodeType dn = Lucee::Field<NDIM, T>::writeToFile(io, node,
        "UnstructGridField");
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

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::writeToTxtFile(std::ofstream& txtFl)
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
      grid->setIndex(idx[0]);
// get cell-center coordinates
      grid->getCentroid(xc);

// write coordinates of cell-center
      for (unsigned i = 0; i < NDIM; ++i)
        txtFl << xc[i] << " ";
// now write out actual data at this location
      for (unsigned i = 0; i < this->getNumComponents(); ++i)
        txtFl << rPtr[i] << " ";
      txtFl << std::endl;
    }
  }

  template<unsigned NDIM, typename T>
  TxIoNodeType UnstructGridField<NDIM, T>::readFromFile(TxIoBase& io,
      TxIoNodeType& node, const std::string& nm)
  {
    TxIoNodeType dn = Lucee::Field<NDIM, T>::readFromFile(io, node,
        "UnstructGridField");
    return dn;
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::sync()
  {
#if 0

// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();
    int rank = comm->getRank();// current rank
    int lg[NDIM], ug[NDIM];
    this->fillWithGhosts(lg, ug);// ghost cells

// initiate (non-blocking) receive of data from our neighbors
    startRecv(rank, lg, ug);
// send data over
    send(rank, lg, ug);
// finish receiving data as recieve was non-blocking
    finishRecv();
#endif
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::startRecv(unsigned rank, int lg[NDIM],
      int ug[NDIM])
  {
#if 0

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
      int tag = 1000;
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
      continue;// nothing do if no intersection

      TxMsgStatus ms = comm->template
      startRecv<T>(this->getNumComponents()*recvRgn.getVolume(), srcrank, tag);
      msgStatus[*niter] = ms;// key in msgStatus map is region we expect receive
    }
    isReceiving = true;
#endif
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::send(unsigned rank, int lg[NDIM],
      int ug[NDIM])
  {
#if 0
// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();

    Lucee::Region<NDIM, int> localBox = this->getRegion();

    std::vector<unsigned> neigh = grid->getSendNeighbors(rank, lg, ug);
    std::vector<unsigned>::const_iterator niter;
    for (niter = neigh.begin(); niter != neigh.end(); ++niter)
    {
      int destrank = grid->getRgnRank(*niter);
      Lucee::Region<NDIM, int> otherBox = grid->getNeighborRgn(*niter);
      int tag = 1000;

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
        Lucee::Region<NDIM, int> sendRgn = localBox.intersect(fakeBox.extend(lg, ug));
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
#endif
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::finishRecv()
  {
#if 0
    if (!isReceiving) return; // nothing to receive

// get hold of comm pointer to do all parallel messaging
    TxCommBase *comm = this->getComm();

// complete receives and copy data into ghost cells
    std::map<int, TxMsgStatus>::iterator i = msgStatus.begin();
    for (; i != msgStatus.end(); ++i)
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
#endif
  }

  template<unsigned NDIM, typename T>
  void UnstructGridField<NDIM, T>::appendLuaCallableMethods(
      Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    Field<NDIM, T>::appendLuaCallableMethods(lfm);
// now append local methods
    lfm.appendFunc("set", luaSet);
    lfm.appendFunc("alias", luaAlias);
    lfm.appendFunc("duplicate", luaDuplicate);
    lfm.appendFunc("applyFuncBc", luaSetGhost);
  }

  template<unsigned NDIM, typename T>
  int UnstructGridField<NDIM, T>::luaSet(lua_State *L)
  {
    UnstructGridField<NDIM, T> *sgf = Lucee::PointerHolder<
        UnstructGridField<NDIM, T> >::getObjAsDerived(L);

    if (!lua_isfunction(L, 2))
    {
      Lucee::Except lce(
          "UnstructGridField::luaSet: Must provide a Lua function to 'set' method");
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_pop(L, 1); // WHY IS THIS NEEDED? THIS LOOKS LIKE A POTENTIAL STACK CORRUPTION

    sgf->setFromLuaFunction(L, fnRef);
    return 0;
  }

  template<unsigned NDIM, typename T>
  int UnstructGridField<NDIM, T>::luaSetGhost(lua_State *L)
  {
    UnstructGridField<NDIM, T> *sgf = Lucee::PointerHolder<
        UnstructGridField<NDIM, T> >::getObjAsDerived(L);

    if (!lua_isnumber(L, 3))
    {
      Lucee::Except lce(
          "UnstructGridField::luaSetGhost: Must provide a number to 'applyFuncBc' method");
      throw lce;
    }

// determine direction in which to apply function BCs
    int dir = (int) lua_tonumber(L, 3);
    if (dir < 0 || dir >= NDIM)
    { // incorrect direction specified
      Lucee::Except lce(
          "UnstructGridField::luaSetGhost: Direction must be one of ");
      for (unsigned i = 0; i < NDIM - 1; ++i)
        lce << i << ", ";
      lce << NDIM - 1 << ".";
      lce << " '" << dir << "' specified instead" << std::endl;
      throw lce;
    }

    if (!lua_isstring(L, 4))
    {
      Lucee::Except lce(
          "UnstructGridField::luaFuncBc: Must provide a side to 'applyFuncBc' method.");
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
      throw Lucee::Except(
          "UnstructGridField::luaFuncBc: side should be one of \"lower\" or \"upper\".");

    if (!lua_isfunction(L, 2))
    {
      Lucee::Except lce(
          "UnstructGridField::luaSetGhost: Must provide a Lua function to 'applyFuncBc' method");
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);

    sgf->setGhostFromLuaFunction(L, fnRef, dir, side);
    return 0;
  }

  template<unsigned NDIM, typename T>
  int UnstructGridField<NDIM, T>::luaAlias(lua_State *L)
  {
#if 0
    UnstructGridField<NDIM, T> *sgf
    = Lucee::PointerHolder<UnstructGridField<NDIM, T> >::getObjAsDerived(L);
    if (! lua_isnumber(L, 2) && ! lua_isnumber(L, 3))
    {
      Lucee::Except lce("UnstructGridField::luaAlias: 'alias' method must be passed half-open interval as [s,e)");
      throw lce;
    }
    unsigned s = lua_tonumber(L, 2);
    unsigned e = lua_tonumber(L, 3);

    Lucee::Field<NDIM, T> aliasFld = sgf->getSubCompView(s, e);
    Lucee::UnstructGridField<NDIM, T>* aliasSgf =
    new Lucee::UnstructGridField<NDIM, T>(aliasFld, *sgf->grid);
// copy data for ghost I/O
    for (unsigned d=0; d<NDIM; ++d)
    {
      aliasSgf->lowerWriteGhost[d] = sgf->lowerWriteGhost[d];
      aliasSgf->upperWriteGhost[d] = sgf->upperWriteGhost[d];
    }

    size_t nbytes = sizeof(Lucee::PointerHolder<UnstructGridField<NDIM, T> >);
    Lucee::PointerHolder<UnstructGridField<NDIM, T> > *ph =
    (Lucee::PointerHolder<UnstructGridField<NDIM, T> >*) lua_newuserdata(L, nbytes);
    ph->pointer = aliasSgf;
    ph->pointer->template setBaseType<Lucee::DataStructIfc>();
    ph->pointer->template setDerivedType<UnstructGridField<NDIM, T> >();

    luaL_getmetatable(L, typeid(UnstructGridField<NDIM, T>).name());
    lua_setmetatable(L, -2);
#endif
    return 1;
  }


  template<unsigned NDIM, typename T>
  int UnstructGridField<NDIM, T>::luaDuplicate(lua_State *L)
  {
#if 0
    UnstructGridField<NDIM, T> *sgf = Lucee::PointerHolder<
        UnstructGridField<NDIM, T> >::getObjAsDerived(L);

    Lucee::UnstructGridField<NDIM, T>* aliasSgf = new Lucee::UnstructGridField<
        NDIM, T>(sgf->duplicate(), *sgf->grid);

    aliasSgf->dataLoc = sgf->dataLoc;
    for (unsigned d = 0; d < NDIM; ++d)
    {
      aliasSgf->lowerWriteGhost[d] = sgf->lowerWriteGhost[d];
      aliasSgf->upperWriteGhost[d] = sgf->upperWriteGhost[d];
    }

    size_t nbytes = sizeof(Lucee::PointerHolder<UnstructGridField<NDIM, T> >);
    Lucee::PointerHolder<UnstructGridField<NDIM, T> > *ph =
        (Lucee::PointerHolder<UnstructGridField<NDIM, T> >*) lua_newuserdata(L,
            nbytes);
    ph->pointer = aliasSgf;
    ph->pointer->template setBaseType<Lucee::DataStructIfc>();
    ph->pointer->template setDerivedType<UnstructGridField<NDIM, T> >();

    luaL_getmetatable(L, typeid(UnstructGridField<NDIM, T> ).name());
    lua_setmetatable(L, -2);
#endif
    return 1;
  }

  template <unsigned NDIM, typename T>
    void
    UnstructGridField<NDIM, T>::setFromLuaFunction(lua_State *L, int ref)
    {
#if 0
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
  // push variables on stack
        for (unsigned i=0; i<3; ++i)
          lua_pushnumber(L, xc[i]);
        if (lua_pcall(L, 3, numOut, 0) != 0)
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
#endif
    }

    template <unsigned NDIM, typename T>
    void
    UnstructGridField<NDIM, T>::setGhostFromLuaFunction(lua_State *L, int ref,
      unsigned dir, unsigned side)
    {
#if 0
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
#endif
    }


  template<unsigned NDIM, typename T>
  UnstructGridField<NDIM, T>::UnstructGridField(const Field<NDIM, T>& fld,
      Lucee::UnstructuredGrid<NDIM>& grd) :
      Lucee::Field<NDIM, T>(fld), grid(&grd)
  {
    this->setName(UnstructGridField<NDIM, double>::id);
  }

// instantiations
  template class UnstructGridField<1, int> ;
  template class UnstructGridField<2, int> ;
  template class UnstructGridField<3, int> ;

  template class UnstructGridField<1, float> ;
  template class UnstructGridField<2, float> ;
  template class UnstructGridField<3, float> ;

  template class UnstructGridField<1, double> ;
  template class UnstructGridField<2, double> ;
  template class UnstructGridField<3, double> ;

}
