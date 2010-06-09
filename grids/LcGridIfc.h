/**
 * @file	LcGridIfc.h
 *
 * @brief	Base class for all grids in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_GRID_IFC_H
#define LC_GRID_IFC_H

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
 * Base class for all grids in Lucee.
 */
  class GridIfc : public Lucee::BasicObj
  {
    public:
/** Destroy grid */
      virtual ~GridIfc();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

/**
 * Write grid to given node in HDF5 file.
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
 * Create a new grid. Name should be provided by a derived class.
 */
      GridIfc();
  };
}

#endif // LC_GRID_IFC_H
