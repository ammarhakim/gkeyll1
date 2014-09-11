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
      .append<Lucee::PoissonBracketGyroEquation4D >();
  }
}
