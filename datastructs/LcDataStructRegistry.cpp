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
#include <LcField.h>

namespace Lucee
{
  void
  registerDataStructObjects(Lucee::LuaState& L)
  {
// register grids
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<1, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<2, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<3, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<4, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<5, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<6, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::Field<7, double> >;
  }
}
