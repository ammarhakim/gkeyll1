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
#include <LcLinCombiner.h>
#include <LcMaxwellTm2DUpdater.h>
#include <LcRectCurlUpdater.h>
#include <LcSolverRegistry.h>
#include <LcTXYZFieldSetter.h>
#include <LcWavePropagationUpdater.h>

namespace Lucee
{
  void
  registerSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::MaxwellTm2DUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectCurlUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectCurlUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectCurlUpdater<3> >;
  }
}
