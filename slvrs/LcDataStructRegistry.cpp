/**
 * @file	LcDataStructRegistry.cpp
 *
 * @brief	Class for registering data-struct objects.
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
#include <LcDataStructRegistry.h>
#include <LcFieldFactory.h>

namespace Lucee
{
  void
  registerDataStructObjects(Lucee::LuaState& L)
  {
// register grids
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::DataStructIfc>, 
      Lucee::FieldFactory<1> >;
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::DataStructIfc>, 
      Lucee::FieldFactory<2> >;
    new Lucee::ObjRegistry<Lucee::GenericFactory<Lucee::DataStructIfc>, 
      Lucee::FieldFactory<3> >;
  }
}
