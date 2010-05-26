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
#include <LcObjCreator.h>
#include <LcRegisterModules.h>
#include <LcRteRegistry.h>
#include <LcSolverIfc.h>
#include <LcSolverRegistry.h>

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
    Lucee::registerRteObjects(L);

// register modules into Lua
    Lucee::ObjCreator<Lucee::SolverIfc>::registerModule(L);
  }
}
