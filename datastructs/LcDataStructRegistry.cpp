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
#include <LcStructGridField.h>

namespace Lucee
{
  void
  registerDataStructObjects(Lucee::LuaState& L)
  {
// register grids
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<1, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<2, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<3, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<4, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<5, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<6, double> >;
    new Lucee::ObjRegistry<Lucee::DataStructIfc, Lucee::StructGridField<7, double> >;
  }
}
