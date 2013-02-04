/**
 * @file	LcLibRegistry.cpp
 *
 * @brief	Method for registering basic library object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcGaussianQuadRule.h>
#include <LcLuaModuleRegistry.h>
#include <LcObjRegistry.h>
#include <LcQuadratureRule.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void
  registerLibObjects(Lucee::LuaState& L)
  {
// register objects
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::QuadratureRule> >
      ::Instance()
      .append<Lucee::GaussianQuadRule>();

// register modules
    Lucee::LuaModuleRegistry<Lucee::QuadratureRule>::registerModule(L);
  }
}
