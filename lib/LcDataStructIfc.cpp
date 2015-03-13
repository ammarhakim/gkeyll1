/**
 * @file	LcDataStructIfc.cpp
 *
 * @brief	Base class for all data in Lucee.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>
#include <LcGlobals.h>
#include <LcPointerHolder.h>

// txbase includes
#include <TxIoBase.h>
#ifdef HAVE_HDF5
#include <TxHdf5Base.h>
#endif

// loki includes
#include <loki/Singleton.h>

// std includes
#include <fstream>
#include <iostream>

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
    BasicObj::readInput(tbl);
  }

  void
  DataStructIfc::write(const std::string& nm, double tCurr)
  {
// output prefix
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string outNm = outPrefix + "_" + nm;
    TxCommBase& comm = *this->getComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->createFile(outNm);
    TxIoNodeType rootGrp = io->openGroup(fn, "/");
    TxIoNodeType td = io->createGroup(rootGrp, "timeData");
    io->writeAttribute(td, "vsType", "time");
    io->writeAttribute(td, "vsStep", 0);
    io->writeAttribute(td, "vsTime", tCurr);

    this->writeToFile(*io, fn, this->getName());
    delete io;
  }

  void
  DataStructIfc::read(const std::string& nm)
  {
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string inNm = outPrefix + "_" + nm;
    TxCommBase& comm = *this->getComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->openFile(inNm, "r");
    TxIoNodeType rootGrp = io->openGroup(fn, "/");
    this->readFromFile(*io, rootGrp, this->getName());
    delete io;
  }

  TxIoNodeType
  DataStructIfc::readFromFile(TxIoBase& io, TxIoNodeType& node, const std::string& nm)
  {
// throw exception if we come here
    throw Lucee::Except("DataStructIfc::readFromFile: Not implemented!");
  }

  void
  DataStructIfc::writeToTxtFile(std::ofstream& txtFl)
  {
// throw exception if we come here
    throw Lucee::Except("DataStructIfc::writeToTxtFile: Not implemented!");
  }

  void
  DataStructIfc::sync()
  {
// throw exception if we come here
    throw Lucee::Except("DataStructIfc::sync: Not implemented!");
  }

  void
  DataStructIfc::appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
  {
    lfm.appendFunc("write", luaWrite);
    lfm.appendFunc("read", luaRead);
    lfm.appendFunc("sync", luaSync);
  }

  int
  DataStructIfc::luaWrite(lua_State *L)
  {
    DataStructIfc *d = Lucee::PointerHolder<DataStructIfc>::getObj(L);
    std::string nm = lua_tostring(L, 2);

    double tCurr = 0.0;
// check if time was specified
    unsigned nArgs = lua_gettop(L);
    if (nArgs == 3)
    {
      if (!lua_isnumber(L, 3))
      {
        Lucee::Except lce("DataStructIfc::luaWrite: Must provide a number when specifying time");
        throw lce;
      }
      tCurr = lua_tonumber(L, 3);
    }
    if (d->isSafeToWrite())
      d->write(nm, tCurr);

    return 0;
  }

  int
  DataStructIfc::luaRead(lua_State *L)
  {
    DataStructIfc *d = Lucee::PointerHolder<DataStructIfc>::getObj(L);
    std::string nm = lua_tostring(L, 2);
    d->read(nm);

    return 0;
  }

  int
  DataStructIfc::luaSync(lua_State *L)
  {
    DataStructIfc *d = Lucee::PointerHolder<DataStructIfc>::getObj(L);
    d->sync();
    return 0;
  }

  DataStructIfc*
  DataStructIfc::clone() const
  {
    return 0;
  }
}
