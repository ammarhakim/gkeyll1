/**
 * @file	LcHyperEquationRegistry.cpp
 *
 * @brief	Method for registering hyperbolic equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcAdvectionEquation.h>
#include <LcEulerEquation.h>
#include <LcHyperEquationRegistry.h>
#include <LcMaxwellEquation.h>

namespace Lucee
{
  void
  registerHyperEquationsObjects(Lucee::LuaState& L)
  {
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AdvectionEquation>;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::EulerEquation>;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::MaxwellEquation>;
  }
}
