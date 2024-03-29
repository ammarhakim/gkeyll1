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
#include <LcNonUniRectCartGrid.h>
#include <LcRectCartGrid.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void
  registerGridObjects(Lucee::LuaState& L)
  {
// register grids
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::GridIfc> >
      ::Instance()
      .append<Lucee::RectCartGrid<1> >()
      .append<Lucee::RectCartGrid<2> >()
      .append<Lucee::RectCartGrid<3> >()
      .append<Lucee::RectCartGrid<4> >()
      .append<Lucee::RectCartGrid<5> >()

      .append<Lucee::NonUniRectCartGrid<1> >()
      .append<Lucee::NonUniRectCartGrid<2> >()
      .append<Lucee::NonUniRectCartGrid<3> >()
      .append<Lucee::NonUniRectCartGrid<4> >()
      .append<Lucee::NonUniRectCartGrid<5> >()

      .append<Lucee::MappedCartGrid<1> >()
      .append<Lucee::MappedCartGrid<2> >()
      .append<Lucee::MappedCartGrid<3> >();
  }
}
