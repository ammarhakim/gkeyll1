/**
 * @file	LcRteRegistry.cpp
 *
 * @brief	Class for registering RTE solver object.
 */

// lucee includes
#include <LcHgPhaseFunction.h>
#include <LcLuaModuleRegistry.h>
#include <LcPlCoeffsPhaseFunction.h>
#include <LcRegisteredObjList.h>
#include <LcRteHomogeneousSlab.h>
#include <LcRtePhaseFunction.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void registerRteObjects(Lucee::LuaState& L)
  {
// register RTE solvers
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::SolverIfc> >
      ::Instance()
      .append<Lucee::RteHomogeneousSlab>();

// register phase functions
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::RtePhaseFunction> >
      ::Instance()
      .append<Lucee::HgPhaseFunction>()
      .append<Lucee::PlCoeffsPhaseFunction>();

// register phase function library into Lua
    Lucee::LuaModuleRegistry<Lucee::RtePhaseFunction>::registerModule(L);
  }
}

