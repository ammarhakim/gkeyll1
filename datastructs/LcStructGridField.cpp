/**
 * @file	LcStructGridField.cpp
 *
 * @brief	StructGridFields are fields that live on structured grids.
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
#include <LcFieldFactory.h>
#include <LcPointerHolder.h>
#include <LcStructGridField.h>

// std includes
#include <fstream>

namespace Lucee
{
// names used in registration system
  template <> const char *StructGridField<1, double>::id = "Field1D";
  template <> const char *StructGridField<2, double>::id = "Field2D";
  template <> const char *StructGridField<3, double>::id = "Field3D";
  template <> const char *StructGridField<4, double>::id = "Field4D";
  template <> const char *StructGridField<5, double>::id = "Field5D";
  template <> const char *StructGridField<6, double>::id = "Field6D";
  template <> const char *StructGridField<7, double>::id = "Field7D";

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField()
    : Lucee::Field<NDIM, T>(), grid(0)
  {
    this->setName(StructGridField<NDIM, double>::id);
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(Lucee::StructuredGridBase<NDIM>* grid, unsigned nc,
        int lg[NDIM], int ug[NDIM])
    : Lucee::Field<NDIM, T>(grid->getLocalBox(), nc, lg, ug, (T) 0),
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
  StructGridField<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::FieldFactory<NDIM, T> ff;
    ff.readInput(tbl);
// re-rest ourself from factory produced grid
    Field<NDIM, T>* nf = ff.create();
    Field<NDIM, T>::operator=(*nf);
    grid = ff.getGridPtr(); // set our grid pointer
    delete nf; // is this really needed?
  }

  template <unsigned NDIM, typename T>
  Lucee::IoNodeType
  StructGridField<NDIM, T>::writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node, const std::string& nm)
  {
// first write the field data to file
    Lucee::IoNodeType dn 
      = Lucee::Field<NDIM, T>::writeToFile(io, node, "StructGridField");
// annotate with viz-schema marks
    io.writeStrAttribute(dn, "vsType", "variable");
    io.writeStrAttribute(dn, "vsMesh", "StructGrid");
    io.writeStrAttribute(dn, "vsCentering", "zonal");
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
      grid->getCentriod(xc);

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
    lua_pop(L, 1);

    sgf->setFromLuaFunction(L, fnRef);
    return 0;
  }

  template <unsigned NDIM, typename T>
  int
  StructGridField<NDIM, T>::luaSetGhost(lua_State *L)
  {
    StructGridField<NDIM, T> *sgf
      = Lucee::PointerHolder<StructGridField<NDIM, T> >::getObjAsDerived(L);
    if (! lua_isfunction(L, 2))
    {
      Lucee::Except lce(
        "StructGridField::luaSetGhost: Must provide a Lua function to 'applyFuncBc' method");
      throw lce;
    }
    int fnRef = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_pop(L, 1);

    if (! lua_isnumber(L, 3))
    {
      Lucee::Except lce(
        "StructGridField::luaSetGhost: Must provide a number to 'applyFuncBc' method");
      throw lce;
    }
// determine direction in which to apply copy BCs
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
      grid->getCentriod(xc); // cell center coordinate
// push function object on stack
      lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
// push variables on stack
      for (unsigned i=0; i<3; ++i)
        lua_pushnumber(L, xc[i]);
      if (lua_pcall(L, 3, numOut, 0) != 0)
      {
        Lucee::Except lce("StructGridField::setFromLuaFunction: ");
        lce << "Problem evaluating function supplied to 'set' method";
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
      up[i] = this->getGlobalUpperExt (i);
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
      grid->getCentriod(xc); // cell center coordinate
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
  template class StructGridField<6, int>;
  template class StructGridField<7, int>;

  template class StructGridField<1, float>;
  template class StructGridField<2, float>;
  template class StructGridField<3, float>;
  template class StructGridField<4, float>;
  template class StructGridField<5, float>;
  template class StructGridField<6, float>;
  template class StructGridField<7, float>;

  template class StructGridField<1, double>;
  template class StructGridField<2, double>;
  template class StructGridField<3, double>;
  template class StructGridField<4, double>;
  template class StructGridField<5, double>;
  template class StructGridField<6, double>;
  template class StructGridField<7, double>;
}
