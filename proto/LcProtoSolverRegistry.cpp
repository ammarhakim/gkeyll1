/**
 * @file	LcProtoSolverRegistry.cpp
 *
 * @brief	Method for registering basic solver object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMaxwellTm2DUpdater.h>
#include <LcMusclHancock1DUpdater.h>
#include <LcProtoSolverRegistry.h>

namespace Lucee
{
  void
  registerProtoSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::MaxwellTm2DUpdater>;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::MusclHancock1DUpdater>;
  }
}
