/**
 * @file	LcSolverRegistry.cpp
 *
 * @brief	Method for registering basic solver object.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcSolverRegistry.h>
#include <LcSolverAssembly.h>

namespace Lucee
{
  void
  registerSolverObjects(Lucee::LuaState& L)
  {
// register solver assembly
    new Lucee::ObjRegistry<Lucee::SolverIfc, Lucee::SolverAssembly>;
  }
}
