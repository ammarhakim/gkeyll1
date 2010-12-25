/**
 * @file	LcHyperEquationRegistry.cpp
 *
 * @brief	Method for registering hyperbolic equations.
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
#include <LcAdvectionEquation.h>
#include <LcEulerEquation.h>
#include <LcHyperEquationRegistry.h>

namespace Lucee
{
  void
  registerHyperEquationsObjects(Lucee::LuaState& L)
  {
// register Euler equations
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::EulerEquation>;
// register advection equation
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AdvectionEquation>;
  }
}
