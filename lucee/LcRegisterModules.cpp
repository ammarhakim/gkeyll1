/**
 * @file	LcRegisterModules.cpp
 *
 * @brief	Class for registering all modules.
 *
 * @version	$Id: LcRegisterModules.cpp 321 2010-03-04 23:13:16Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// lucee includes
#include <LcDataStructIfc.h>
#include <LcDataStructRegistry.h>
#include <LcFunctionIfc.h>
#include <LcGenericFactory.h>
#include <LcGridIfc.h>
#include <LcGridRegistry.h>
#include <LcHyperEquation.h>
#include <LcHyperEquationRegistry.h>
#include <LcLibRegistry.h>
#include <LcLuaModuleRegistry.h>
#include <LcRegisterModules.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>
#include <LcSolverRegistry.h>
#include <LcUpdaterIfc.h>

namespace Lucee
{
/**
 * Register all modules in Lucee.
 *
 * @param L Lua state object in which modules should be registered.
 */
  void registerModules(Lucee::LuaState& L)
  {
// register objects
    Lucee::registerSolverObjects(L);
    Lucee::registerHyperEquationsObjects(L);
    Lucee::registerRteObjects(L);
    Lucee::registerGridObjects(L);
    Lucee::registerDataStructObjects(L);
    Lucee::registerLibObjects(L);

// register modules into Lua
    Lucee::LuaModuleRegistry<Lucee::SolverIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::HyperEquation>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::GridIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::DataStructIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::UpdaterIfc>::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::FunctionIfc>::registerModule(L);
  }
}
