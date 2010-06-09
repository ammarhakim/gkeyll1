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
#include <LcGridRegistry.h>
#include <LcRectCartGridFactory.h>

namespace Lucee
{
  void
  registerGridObjects(Lucee::LuaState& L)
  {
// register grids
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::GridIfc>, 
      Lucee::RectCartGridFactory<1> >;
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::GridIfc>, 
      Lucee::RectCartGridFactory<2> >;
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::GridIfc>, 
      Lucee::RectCartGridFactory<3> >;
  }
}
