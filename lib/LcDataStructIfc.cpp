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
// delete I/O pointer
      delete io;
    }
    else
    {
// open text file to write in
      std::ofstream outFl(outNm.c_str(), std::fstream::out);
// write out data
      this->writeToTxtFile(outFl);
    }
  }

  void
  DataStructIfc::writeToTxtFile(std::ofstream& txtFl)
  {
// throw exception if we come here
    throw Lucee::Except("DataStructIfc::writeToTxtFile: Not implemented!");
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
