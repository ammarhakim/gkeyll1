/**
 * @file	LcDataStructIfc.h
 *
 * @brief	Base class for all data in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_DATA_STRUCT_IFC_H
#define LC_DATA_STRUCT_IFC_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcIoBase.h>
#include <LcLuaTable.h>

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
/** Destroy dataStruct */
      virtual ~DataStructIfc();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Clone data and return pointer to cloned object.
 *
 * @return cloned object.
 */
      virtual DataStructIfc* clone() const;

/**
 * Write dataStruct to given node in HDF5 file.
 *
 * @param io I/O object for I/O.
 * @param node Node to write to.
 * @param nm Name of the grid as it should appear in output.
 * @return node to which data was written.
 */
      virtual Lucee::IoNodeType writeToFile(Lucee::IoBase& io, Lucee::IoNodeType& node,
        const std::string& nm) = 0;

    protected:
/**
 * Create a new dataStruct. Name should be provided by a derived
 * class.
 */
      DataStructIfc();
  };
}

#endif // LC_DATA_STRUCT_IFC_H
