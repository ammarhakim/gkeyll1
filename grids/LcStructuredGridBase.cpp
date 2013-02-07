/**
 * @file	LcStructuredGridBase.cpp
 *
 * @brief	Base class for body fitted grid in arbitrary dimensions.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcPointerHolder.h>
#include <LcStructuredGridBase.h>

// txbase includes
#include <TxCommBase.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  template <unsigned NDIM>
  StructuredGridBase<NDIM>::~StructuredGridBase()
  {
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::readInput(Lucee::LuaTable& tbl)
  {
// read number of cells in global region
    std::vector<double> cells = tbl.getNumVec("cells");
    if (cells.size() != NDIM)
    {
      Lucee::Except lce("StructuredGridBase::readInput: 'cells' should have exactly");
      lce << NDIM << " elements. Instead has " << cells.size() << std::endl;
      throw lce;
    }

// create global region
    int zeros[NDIM], icells[NDIM];
    for (unsigned i=0; i<NDIM; ++i)
    {
      zeros[i] = 0;
      icells[i] = (int) cells[i];
    }
    globalRgn = Lucee::Region<NDIM, int>(zeros, icells);

    std::vector<double> lower(NDIM), upper(NDIM);
// read lower limits of computational space
    if (tbl.hasNumVec("lower"))
    {
      lower = tbl.getNumVec("lower");
      if (lower.size() != NDIM)
      {
        Lucee::Except lce("StructuredGridBase::readInput: 'lower' should have exactly ");
        lce << NDIM << " elements. Instead has " << lower.size() << std::endl;
        throw lce;
      }
    }
    else
    {
      for (unsigned i=0; i<NDIM; ++i) lower[i] = 0.0;
    }

// read upper limits of computational space
    if (tbl.hasNumVec("upper"))
    {
      upper = tbl.getNumVec("upper");
      if (upper.size() != NDIM)
      {
        Lucee::Except lce("StructuredGridBase::readInput: 'upper' should have exactly ");
        lce << NDIM << " elements. Instead has " << upper.size() << std::endl;
        throw lce;
      }
    }
    else
    {
      for (unsigned i=0; i<NDIM; ++i) upper[i] = 1.0;
    }
    compSpace = Lucee::Region<NDIM, double>(&lower[0], &upper[0]);

// get comm pointer
    TxCommBase *comm = Loki::SingletonHolder<Lucee::Globals>
      ::Instance().comm;

// create decomposed region
    decompRgn.reset( new Lucee::DecompRegion<NDIM>(globalRgn) );

#ifdef HAVE_MPI
// get decomposition object
    if (tbl.template hasObject<Lucee::DecompRegionCalcIfc<NDIM> >("decomposition"))
    {
      Lucee::DecompRegionCalcIfc<NDIM>& decompCalc
        = tbl.template getObject<Lucee::DecompRegionCalcIfc<NDIM> >("decomposition");
// compute decomposition
      decompCalc.calcDecomp(comm->getNumProcs(), *decompRgn);
    }
#endif
// compute local region
    localRgn = decompRgn->getRegion(comm->getRank());
  }

  template <unsigned NDIM>
  unsigned
  StructuredGridBase<NDIM>::getNumCells(unsigned dir) const
  {
    return globalRgn.getShape(dir);
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getGlobalRegion() const 
  { 
    return globalRgn; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, int>
  StructuredGridBase<NDIM>::getLocalRegion() const 
  { 
    return localRgn; 
  }

  template <unsigned NDIM>
  Lucee::Region<NDIM, double>
  StructuredGridBase<NDIM>::getComputationalSpace() const 
  { 
    return compSpace; 
  }

  template <unsigned NDIM>
  double
  StructuredGridBase<NDIM>::getDx(unsigned dir) const
  { 
    return compSpace.getShape(dir)/globalRgn.getShape(dir);
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
  std::vector<unsigned> 
  StructuredGridBase<NDIM>::getRecvNeighbors(unsigned rn, 
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    return decompRgn->getRecvNeighbors(rn, lowerExt, upperExt);
  }

  template <unsigned NDIM>
  std::vector<unsigned> 
  StructuredGridBase<NDIM>::getSendNeighbors(unsigned rn, 
    const int lowerExt[NDIM], const int upperExt[NDIM]) const
  {
    return decompRgn->getSendNeighbors(rn, lowerExt, upperExt);
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
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> lb = g->getLocalRegion();
    lua_pushnumber(L, lb.getLower(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetLocalUpper(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> lb = g->getLocalRegion();
    lua_pushnumber(L, lb.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetGlobalLower(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> gb = g->getGlobalRegion();
    lua_pushnumber(L, gb.getLower(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetGlobalUpper(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, int> gb = g->getGlobalRegion();
    lua_pushnumber(L, gb.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  int
  StructuredGridBase<NDIM>::luaGetLowerCoord(lua_State *L)
  {
    StructuredGridBase<NDIM> *g
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

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
      = Lucee::PointerHolder<StructuredGridBase<NDIM> >::getObj(L);

    int dir = lua_tonumber(L, 2);
    Lucee::Region<NDIM, double> cs = g->getComputationalSpace();
    lua_pushnumber(L, cs.getUpper(dir));
    return 1;
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase()
    : localRgn(&Lucee::FixedVector<NDIM, int>(1)[0]),
      globalRgn(&Lucee::FixedVector<NDIM, int>(1)[0]),
      compSpace(&Lucee::FixedVector<NDIM, double>(1.0)[0]),
      decompRgn( new Lucee::DecompRegion<NDIM>(globalRgn) )
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase(const Lucee::Region<NDIM, int>& globalRgn,
    const Lucee::Region<NDIM, double>& compSpace)
    : localRgn(globalRgn), globalRgn(globalRgn), compSpace(compSpace),
      decompRgn( new Lucee::DecompRegion<NDIM>(globalRgn) )
  {
  }

  template <unsigned NDIM>
  StructuredGridBase<NDIM>::StructuredGridBase(const Lucee::DecompRegion<NDIM>& dcmpRgn,
    const Lucee::Region<NDIM, double>& compSpace)
    : localRgn(dcmpRgn.getRegion(
        Loki::SingletonHolder<Lucee::Globals>
        ::Instance().comm->getRank()
                                 )
               ),
      globalRgn(dcmpRgn.getGlobalRegion()), compSpace(compSpace),
      decompRgn( new Lucee::DecompRegion<NDIM>(dcmpRgn) )
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
    localRgn = sg.localRgn;
    globalRgn = sg.globalRgn;
    compSpace = sg.compSpace;
    decompRgn = sg.decompRgn;
    
    return *this;
  }

  template <unsigned NDIM>
  void
  StructuredGridBase<NDIM>::setGridData(const Lucee::Region<NDIM, int>& lb,
    const Lucee::Region<NDIM, int>& gb, const Lucee::Region<NDIM, double>& cs)
  {
    localRgn = lb;
    globalRgn = gb;
    compSpace = cs;
    decompRgn.reset( new Lucee::DecompRegion<NDIM>(globalRgn) );
  }

// instantiations
  template class StructuredGridBase<1>;
  template class StructuredGridBase<2>;
  template class StructuredGridBase<3>;
  template class StructuredGridBase<4>;
  template class StructuredGridBase<5>;
}
