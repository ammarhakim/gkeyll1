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
  }

  template <unsigned NDIM, typename T>
  StructGridField<NDIM, T>::StructGridField(const StructGridField<NDIM, T>& fld)
    : Lucee::Field<NDIM, T>(fld), grid(fld.grid)
  {
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
  void
  StructGridField<NDIM, T>::readInput(Lucee::LuaTable& tbl)
  {
    Lucee::FieldFactory<NDIM, T> ff;
    ff.readInput(tbl);
// re-rest ourself from factory produced grid
    Field<NDIM, T>::operator=(*ff.create());
    grid = ff.getGridPtr(); // set out grid pointer
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
  StructGridField<NDIM, T>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// call base class Lua methods
    Field<NDIM, T>::appendLuaCallableMethods(lfm);
// now append local methods
    lfm.appendFunc("set", luaSet);
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
