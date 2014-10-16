/**
 * @file	LcPoissonBracketEquationRegistry.cpp
 *
 * @brief	Method for registering poisson bracket equations.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcPoissonBracketEquationRegistry.h>
#include <LcPoissonBracketAdvectionEquation4D.h>
#include <LcPoissonBracketCanonical2D.h>
#include <LcPoissonBracketCanonical4D.h>
#include <LcPoissonBracketGyroEquation4D.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void
  registerPoissonBracketEquationObjects(Lucee::LuaState& L)
  {
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::PoissonBracketEquation> >
      ::Instance()
      .append<Lucee::PoissonBracketAdvectionEquation4D >()
      .append<Lucee::PoissonBracketCanonical2D >()
      .append<Lucee::PoissonBracketCanonical4D >()
      .append<Lucee::PoissonBracketGyroEquation4D >();
  }
}