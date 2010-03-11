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
#include <LcLuaState.h>
#include <LcObjCreator.h>
#include <LcRegisterModules.h>
#include <LcSolverIfc.h>

namespace Lucee
{
  void registerModules(Lucee::LuaState& L)
  {
// register objects
    Lucee::registerRteObjects(L);

// register modules into Lua
    Lucee::ObjCreator<Lucee::SolverIfc>::registerModule(L);
  }
}
