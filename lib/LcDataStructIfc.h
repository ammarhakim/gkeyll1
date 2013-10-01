/**
 * @file	LcDataStructIfc.h
 *
 * @brief	Base class for all data in Lucee.
 */

#ifndef LC_DATA_STRUCT_IFC_H
#define LC_DATA_STRUCT_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcLuaTable.h>

// txbase includes
#include <TxIoBase.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Base class for all dataStruct in Lucee.
 */
  class DataStructIfc : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/** Destroy dataStruct */
      virtual ~DataStructIfc();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Write data-structure to specified HDF5 file.
 *
 * @param nm Name of file to write.
 * @param tm Time which data was written.
 */
      void write(const std::string& nm, double tm);

/**
 * Read data-structure from specified HDF5 file.
 *
 * @param nm Name of file to read.
 */
      void read(const std::string& nm);

/**
 * Clone data and return pointer to cloned object.
 *
 * @return cloned object.
 */
      virtual DataStructIfc* clone() const;

/**
 * Get communicator object to perform I/O.
 *
 * @param Reference to a communicator.
 */
      virtual TxCommBase& getDataComm();

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the data-struct as it should appear in output.
 * @return node to which data was written.
 */
      virtual TxIoNodeType writeToFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm) = 0;

/**
 * Read dataStruct from given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to read data from.
 * @param nm Name of the data-struct as it appears in input.
 * @return node from which data was read.
 */
      virtual TxIoNodeType readFromFile(TxIoBase& io, TxIoNodeType& node,
        const std::string& nm);

/**
 * Write dataStruct to specified text file. The data is written as
 * plain text.
 *
 * @param txtFl Text file handle for output.
 */
      virtual void writeToTxtFile(std::ofstream& txtFl);

/**
 * Synchronise data in ghost cells, copying skin data into neighbor
 * ghost cells.
 */
      virtual void sync();

/**
 * Method that performs registration of Lua functions.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Lua callable method for writing out data to HDF5 file.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaWrite(lua_State *L);

/**
 * Lua callable method for reading data from HDF5 file.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaRead(lua_State *L);

/**
 * Lua callable method for synchronizing ghost cell data.
 *
 * @param L Lua state to use.
 * @return number of output parameters.
 */
      static int luaSync(lua_State *L);

    protected:
/**
 * Create a new dataStruct. Name should be provided by a derived
 * class.
 */
      DataStructIfc();
  };
}

#endif // LC_DATA_STRUCT_IFC_H
