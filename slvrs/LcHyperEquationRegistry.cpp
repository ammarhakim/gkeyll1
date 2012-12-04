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

namespace Lucee
{
  void
  registerHyperEquationsObjects(Lucee::LuaState& L)
  {
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AdvectionEquation>;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::EulerEquation>;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::MaxwellEquation>;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::PhMaxwellEquation>;

    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::DivEquation<1> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::DivEquation<2> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::DivEquation<3> >;

    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::GradEquation<1> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::GradEquation<2> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::GradEquation<3> >;

    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AuxAdvectionEquation<1> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AuxAdvectionEquation<2> >;
    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::AuxAdvectionEquation<3> >;

    new Lucee::ObjRegistry<Lucee::HyperEquation, Lucee::LenardBernsteinVParEquation >;
  }
}
