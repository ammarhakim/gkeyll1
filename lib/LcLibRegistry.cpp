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
#include <LcFunctionIfc.h>
#include <LcLuaTXYZFunction.h>
#include <LcObjRegistry.h>

namespace Lucee
{
  void
  registerLibObjects(Lucee::LuaState& L)
  {
// registry functions
    new Lucee::ObjRegistry<Lucee::FunctionIfc, Lucee::LuaTXYZFunction>;
  }
}
