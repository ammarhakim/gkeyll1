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
#include <LcAuxAdvectionEquation.h>
#include <LcDivEquation.h>
#include <LcEulerEquation.h>
#include <LcGradEquation.h>
#include <LcHyperEquationRegistry.h>
#include <LcLenardBernsteinVParEquation.h>
#include <LcMaxwellEquation.h>
#include <LcPhMaxwellEquation.h>
#include <LcRegisteredObjList.h>

// loki includes
#include <loki/Singleton.h>

namespace Lucee
{
  void
  registerHyperEquationsObjects(Lucee::LuaState& L)
  {
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::HyperEquation> >
      ::Instance()
      .append<Lucee::AdvectionEquation>()
      .append<Lucee::EulerEquation>()
      .append<Lucee::MaxwellEquation>()
      .append<Lucee::PhMaxwellEquation>()
      .append<Lucee::DivEquation<1> >()
      .append<Lucee::DivEquation<2> >()
      .append<Lucee::DivEquation<3> >()
      .append<Lucee::GradEquation<1> >()
      .append<Lucee::GradEquation<2> >()
      .append<Lucee::GradEquation<3> >()
      .append<Lucee::AuxAdvectionEquation<1> >()
      .append<Lucee::AuxAdvectionEquation<2> >()
      .append<Lucee::AuxAdvectionEquation<3> >()
      .append<Lucee::LenardBernsteinVParEquation >();
  }
}
