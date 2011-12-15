/**
 * @file	LcGridRegistry.cpp
 *
 * @brief	Class for registering grid objects.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGridRegistry.h>
#include <LcMappedCartGrid.h>
#include <LcRectCartGrid.h>

namespace Lucee
{
  void
  registerGridObjects(Lucee::LuaState& L)
  {
// register grids
    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::RectCartGrid<1> >;
    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::RectCartGrid<2> >;
    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::RectCartGrid<3> >;

    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::MappedCartGrid<1> >;
    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::MappedCartGrid<2> >;
    new Lucee::ObjRegistry<Lucee::GridIfc, Lucee::MappedCartGrid<3> >;
  }
}
