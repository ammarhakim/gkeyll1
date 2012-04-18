/**
 * @file	LcProtoSolverRegistry.cpp
 *
 * @brief	Method for registering basic solver object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcEnergyFromStreamAndVortUpdater.h>
#include <LcEnergyFromStreamFunctionUpdater.h>
#include <LcEnstrophyUpdater.h>
#include <LcMaxwellTm2DUpdater.h>
#include <LcModalDg1DUpdater.h>
#include <LcModalDgLimiter1DUpdater.h>
#include <LcMusclHancock1DUpdater.h>
#include <LcNodalGradientUpdater.h>
#include <LcNodalPoissonBracketUpdater.h>
#include <LcProtoSolverRegistry.h>
#include <LcRectSecondOrderCentralDiffUpdater.h>

#ifdef HAVE_FFTW3
# include <LcPeriodicPoisson2DUpdater.h>
#endif

namespace Lucee
{
  void
  registerProtoSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::MaxwellTm2DUpdater>;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::MusclHancock1DUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectSecondOrderCentralDiffUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectSecondOrderCentralDiffUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::RectSecondOrderCentralDiffUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::ModalDg1DUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::ModalDgLimiter1DUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalPoissonBracketUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EnergyFromStreamFunctionUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EnergyFromStreamAndVortUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EnstrophyUpdater>;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalGradientUpdater>;

#ifdef HAVE_FFTW3
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::PeriodicPoisson2DUpdater>;
#endif

  }
}
