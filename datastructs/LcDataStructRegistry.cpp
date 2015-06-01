/**
 * @file	LcDataStructRegistry.cpp
 *
 * @brief	Class for registering data-struct objects.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructRegistry.h>
#include <LcDynVector.h>
#include <LcRegisteredObjList.h>
#include <LcStructGridField.h>
#include <LcUnstructGridField.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void registerDataStructObjects(Lucee::LuaState& L)
  {
// register data-structures
    Loki::SingletonHolder < Lucee::RegisteredObjList<Lucee::DataStructIfc>
        > ::Instance().append<Lucee::StructGridField<1, double> >().append<
            Lucee::StructGridField<2, double> >().append<
            Lucee::StructGridField<3, double> >().append<
            Lucee::StructGridField<4, double> >().append<
            Lucee::StructGridField<5, double> >().append<
            Lucee::UnstructGridField<1, double> >().append<
            Lucee::UnstructGridField<2, double> >().append<
            Lucee::UnstructGridField<3, double> >().append<
            Lucee::DynVector<double> >();
  }
}
