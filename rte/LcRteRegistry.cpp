/**
 * @file	LcRteRegistry.cpp
 *
 * @brief	Class for registering RTE solver object.
 */

// lucee includes
#include <LcHgPhaseFunction.h>
#include <LcLuaModuleRegistry.h>
#include <LcPlCoeffsPhaseFunction.h>
#include <LcRteHomogeneousSlab.h>
#include <LcRtePhaseFunction.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>

namespace Lucee
{
  void registerRteObjects(Lucee::LuaState& L)
  {
// register RTE solvers
    new Lucee::ObjRegistry<Lucee::SolverIfc, Lucee::RteHomogeneousSlab>;

// register phase functions
    new Lucee::ObjRegistry<Lucee::RtePhaseFunction, Lucee::HgPhaseFunction>;
    new Lucee::ObjRegistry<Lucee::RtePhaseFunction, Lucee::PlCoeffsPhaseFunction>;

// register phase function library into Lua
    Lucee::LuaModuleRegistry<Lucee::RtePhaseFunction>::registerModule(L);
  }
}

