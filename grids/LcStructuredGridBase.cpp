/**
 * @file	LcStructuredGridBase.cpp
 *
 * @brief	Base class for body fitted grid in arbitrary dimensions.
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
#include <LcPointerHolder.h>
#include <LcStructuredGridBase.h>

namespace Lucee
{
  template <unsigned NDIM>
  StructuredGridBase<NDIM>::~StructuredGridBase()
  {
  }

  template <unsigned NDIM>
  unsigned
  StructuredGridBase<NDIM>::getNumCells(unsigned dir) const
  {
    return globalBox.getShape(dir);
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getGlobalBox() const 
  { 
    return globalBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getLocalBox() const 
  { 
    return localBox; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, double>
  StructuredGridBase<NDIM>::getComputationalSpace() const 
  { 
    return compSpace; 
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i) const
  {
    currIdx[0] = i;
  }
  
  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i, int j) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(int i, int j, int k) const
  {
    currIdx[0] = i;
    currIdx[1] = j;
    currIdx[2] = k;
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setIndex(const int idx[NDIM]) const
  {
    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = idx[i];
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("localLowerIndex", luaGetLocalLower);
    lfm.appendFunc("localUpperIndex", luaGetLocalUpper);
    lfm.appendFunc("globalLowerIndex", luaGetGlobalLower);
    lfm.appendFunc("globalUpperIndex", luaGetGlobalUpper);
    lfm.appendFunc("lower", luaGetLowerCoord);
    lfm.appendFunc("upper", luaGetUpperCoord);
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetLocalLower(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> lb = g->getLocalBox();
    lua_pushnumber(L, lb.getLower(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetLocalUpper(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> lb = g->getLocalBox();
    lua_pushnumber(L, lb.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetGlobalLower(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> gb = g->getGlobalBox();
    lua_pushnumber(L, gb.getLower(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetGlobalUpper(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> gb = g->getGlobalBox();
    lua_pushnumber(L, gb.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetLowerCoord(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, double> cs = g->getComputationalSpace();
    lua_pushnumber(L, cs.getLower(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetUpperCoord(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::checkUserType(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, double> cs = g->getComputationalSpace();
    lua_pushnumber(L, cs.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase()
    : localBox(&Lucee::FixedVector<NDIM, int>(1)[0]),
      globalBox(&Lucee::FixedVector<NDIM, int>(1)[0]),
      compSpace(&Lucee::FixedVector<NDIM, double>(1.0)[0])
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase(const Lucee::Region<NDIM, int>& localBox,
    const Lucee::Region<NDIM, int>& globalBox,
    const Lucee::Region<NDIM, double>& compSpace)
    : localBox(localBox), globalBox(globalBox), compSpace(compSpace)
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>&
  StructuredGridBase<NDIM>::operator=(const StructuredGridBase<NDIM>& sg)
  {
    if (&sg == this)
      return *this;

    for (unsigned i=0; i<NDIM; ++i)
      currIdx[i] = sg.currIdx[i];
    localBox = sg.localBox;
    globalBox = sg.globalBox;
    compSpace = sg.compSpace;
    
    return *this;
  }

// instantiations
  template class StructuredGridBase<1>;
  template class StructuredGridBase<2>;
  template class StructuredGridBase<3>;
  template class StructuredGridBase<4>;
  template class StructuredGridBase<5>;
  template class StructuredGridBase<6>;
  template class StructuredGridBase<7>;
}
