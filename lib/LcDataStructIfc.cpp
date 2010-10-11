/**
 * @file	LcDataStructIfc.cpp
 *
 * @brief	Base class for all data in Lucee.
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
#include <LcDataStructIfc.h>
#include <LcHdf5Io.h>
#include <LcPointerHolder.h>

namespace Lucee
{
// set madule name
  const char *DataStructIfc::id = "DataStruct";

  DataStructIfc::DataStructIfc()
    : Lucee::BasicObj("DataStruct")
  {
  }

  DataStructIfc::~DataStructIfc()
  {
  }

  void
  DataStructIfc::readInput(Lucee::LuaTable& tbl)
  {
  }

  void
  DataStructIfc::write(const std::string& nm)
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

  void
  DataStructIfc::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("write", luaWrite);
  }

  int
  DataStructIfc::luaWrite(lua_State *L)
  {
    DataStructIfc *d
      = Lucee::PointerHolder<DataStructIfc>::getObj(L);
    std::string nm = lua_tostring(L, 2);
    d->write(nm);

    return 0;
  }

  DataStructIfc*
  DataStructIfc::clone() const
  {
    return 0;
  }
}
