/**
 * @file	LcSolverRegistry.cpp
 *
 * @brief	Method for registering basic solver object.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBcUpdater.h>
#include <LcCopyBoundaryCondition.h>
#include <LcEdgeFaceCurlUpdater.h>
#include <LcFaceEdgeCurlUpdater.h>
#include <LcLinCombiner.h>
#include <LcLuaModuleRegistry.h>
#include <LcSolverRegistry.h>
#include <LcWavePropagationUpdater.h>
#include <LcZeroNormalBoundaryCondition.h>

namespace Lucee
{
  void
  registerSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::LinCombiner<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::WavePropagationUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FaceEdgeCurlUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FaceEdgeCurlUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FaceEdgeCurlUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EdgeFaceCurlUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EdgeFaceCurlUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EdgeFaceCurlUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::BcUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::BcUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::BcUpdater<3> >;

// register boundary conditions
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::CopyBoundaryCondition>;
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::ZeroNormalBoundaryCondition>;

// register boundary condition library into Lucee (this needs to be
// done once here as boundary conditions are local to the slvr
// directory. Perhaps could have also done in lucee directory)
    Lucee::LuaModuleRegistry<Lucee::BoundaryCondition>::registerModule(L);
  }
}
