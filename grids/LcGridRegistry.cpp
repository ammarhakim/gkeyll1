/**
 * @file	LcGridRegistry.cpp
 *
 * @brief	Class for registering grid objects.
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
#include <LcGridRegistry.h>
#include <LcRectCartGrid.h>

namespace Lucee
{
  void
  registerGridObjects(Lucee::LuaState& L)
  {
// register grids
//     new Lucee::ObjRegistry<Lucee::GridBase, Lucee::RectCartGrid<1> >;
//     new Lucee::ObjRegistry<Lucee::GridBase, Lucee::RectCartGrid<2> >;
//     new Lucee::ObjRegistry<Lucee::GridBase, Lucee::RectCartGrid<3> >;
  }
}
