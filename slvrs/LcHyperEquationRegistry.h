/**
 * @file	LcHyperEquationsRegistry.h
 *
 * @brief	Method for registering hyperbolic equations.
 */

#ifndef LC_HYPER_EQUATION_REGISTRY_H
#define LC_HYPER_EQUATION_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register hyperbolic equations.
 */
  void registerHyperEquationsObjects(Lucee::LuaState& L);
}

#endif // LC_HYPER_EQUATION_REGISTRY_H
