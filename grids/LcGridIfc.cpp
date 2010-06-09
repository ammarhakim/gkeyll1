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

namespace Lucee
{
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
}
