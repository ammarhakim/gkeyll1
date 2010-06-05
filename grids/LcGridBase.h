/**
 * @file	LcGridBase.h
 *
 * @brief	Base class for all grids in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_GRID_BASE_H
#define LC_GRID_BASE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>
#include <LcLuaTable.h>

namespace Lucee
{
/**
 * Base class for all grids in Lucee.
 */
  class GridBase : public Lucee::BasicObj
  {
    public:
/** Destroy grid */
      virtual ~GridBase();

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl);

    protected:
/**
 * Create a new grid. Name should be provided by a derived class.
 */
      GridBase();
  };
}

#endif // LC_GRID_BASE_H
