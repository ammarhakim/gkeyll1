/**
 * @file	LcPoissonBracketEquationRegistry.h
 *
 * @brief	Method for registering poisson bracket equations.
 */

#ifndef LC_POISSON_BRACKET_EQUATION_REGISTRY_H
#define LC_POISSON_BRACKET_EQUATION_REGISTRY_H

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
  void registerPoissonBracketEquationObjects(Lucee::LuaState& L);
}

#endif // LC_POISSON_BRACKET_EQUATION_REGISTRY_H
