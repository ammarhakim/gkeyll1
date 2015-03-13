/**
 * @file	LcGridIfc.cpp
 *
 * @brief	Base class for all grids in Lucee.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGlobals.h>
#include <LcGridIfc.h>
#include <LcPointerHolder.h>

// txbase includes
#include <TxHdf5Base.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
// set ids for grid
  const char *Lucee::GridIfc::id = "Grid";

  GridIfc::GridIfc()
    : Lucee::BasicObj("Grid")
  {
  }

  GridIfc::~GridIfc()
  {
  }

  void
  GridIfc::readInput(Lucee::LuaTable& tbl)
  {
    BasicObj::readInput(tbl);
  }

  void
  GridIfc::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
// function to write grid to HDF5 file
    lfm.appendFunc("write", luaWrite);
  }

  int
  GridIfc::luaWrite(lua_State *L)
  {
    GridIfc *g
      = Lucee::PointerHolder<GridIfc>::getObj(L);
    std::string nm = lua_tostring(L, 2);
    g->write(nm);

    return 0;
  }

  void
  GridIfc::write(const std::string& nm)
  {
// output prefix
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string outNm = outPrefix + "_" + nm;
    TxCommBase& comm = *this->getDataComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->createFile(outNm);
    this->writeToFile(*io, fn, this->getName());
    delete io;
  }
}
