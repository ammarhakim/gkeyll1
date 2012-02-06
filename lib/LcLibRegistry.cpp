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

namespace Lucee
{
  void
  registerLibObjects(Lucee::LuaState& L)
  {
// register objects
    new Lucee::ObjRegistry<Lucee::QuadratureRule, Lucee::GaussianQuadRule>;

// register modules
    Lucee::LuaModuleRegistry<Lucee::QuadratureRule>::registerModule(L);
  }
}
