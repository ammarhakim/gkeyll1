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
    bool isH5 = true;
// output prefix
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
// check file extention to determine if hdf5 or plain text should be written
    std::string snm = nm;
    unsigned trunc = nm.find_last_of(".", snm.size());
    if (trunc > 0)
      snm.erase(0, trunc);
    if (snm == ".txt")
      isH5 = false; // write out text file

// output name
    std::string outNm = outPrefix + "_" + nm;

    if (isH5)
    {
#ifdef HAVE_MPI
      TxMpiBase commBase;
#else
      TxSelfBase commBase;
#endif
      TxIoBase *io = new TxHdf5Base(&commBase);
// open file to write in
      TxIoNodeType fn = io->createFile(outNm);
// write data
      this->writeToFile(*io, fn, this->getName());
    }
    else
    {
      throw Lucee::Except("Lucee::GridIfc: Text output of grid not supported yet.");
    }
  }
}
