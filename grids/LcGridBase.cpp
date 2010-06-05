/**
 * @file	LcGridBase.cpp
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
#include <LcGridBase.h>

namespace Lucee
{
// set madule name
  const char *GridBase::id = "Grid";

  GridBase::GridBase()
    : Lucee::BasicObj("Grid")
  {
  }

  GridBase::~GridBase()
  {
  }

  void
  GridBase::readInput(Lucee::LuaTable& tbl)
  {
  }
}
