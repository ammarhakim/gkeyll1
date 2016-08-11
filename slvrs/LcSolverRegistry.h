/**
 * @file	LcSolverRegistry.h
 *
 * @brief	Method for registering basic solver object.
 */

#ifndef LC_SOLVER_REGISTRY_H
#define LC_SOLVER_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register basic solver objects and modules.
 */
  void registerSolverObjects(Lucee::LuaState& L);
}

#endif // LC_SOLVER_REGISTRY_H
