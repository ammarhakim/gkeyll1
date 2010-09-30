/**
 * @file	LcGridIfc.cpp
 *
 * @brief	Base class for all grids in Lucee.
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
#include <LcGridIfc.h>
#include <LcHdf5Io.h>
#include <LcPointerHolder.h>

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
      = Lucee::PointerHolder<GridIfc>::checkUserType(L);
    std::string nm = lua_tostring(L, 2);
    g->write(nm);

    return 0;
  }

  void
  GridIfc::write(const std::string& nm)
  {
#ifdef HAVE_MPI
    Lucee::Hdf5Io io(MPI_COMM_WORLD, MPI_INFO_NULL);
#else
    Lucee::Hdf5Io io(0, 0);
#endif
// open file to write in
    Lucee::IoNodeType fn = io.createFile(nm);
// write data
    this->writeToFile(io, fn, this->getName());
  }
}
