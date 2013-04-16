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
  }

  void
  DataStructIfc::write(const std::string& nm)
  {
// output prefix
    std::string outPrefix = Loki::SingletonHolder<Lucee::Globals>::Instance().outPrefix;
    std::string outNm = outPrefix + "_" + nm;
    TxCommBase& comm = this->getDataComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->createFile(outNm);
    this->writeToFile(*io, fn, this->getName());
    delete io;
  }

  void
  DataStructIfc::read(const std::string& nm, const std::string& grp)
  {
    TxCommBase& comm = this->getDataComm();

    TxIoBase *io = new TxHdf5Base(&comm);
    TxIoNodeType fn = io->openFile(nm, "r");
    this->readFromFile(*io, fn, nm);
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
    lfm.appendFunc("sync", luaSync);
  }

  int
  DataStructIfc::luaWrite(lua_State *L)
  {
    DataStructIfc *d = Lucee::PointerHolder<DataStructIfc>::getObj(L);
    std::string nm = lua_tostring(L, 2);
    d->write(nm);
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

  TxCommBase&
  DataStructIfc::getDataComm()
  {
    return *Loki::SingletonHolder<Lucee::Globals>::Instance().comm;
  }
}
