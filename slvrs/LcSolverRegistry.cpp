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
#include <LcConstBoundaryCondition.h>
#include <LcCopyBoundaryCondition.h>
#include <LcCopyContToDisContFieldUpdater.h>
#include <LcCurrentSource.h>
#include <LcEdgeFaceCurlUpdater.h>
#include <LcEvalOnNodesUpdater.h>
#include <LcFaceEdgeCurlUpdater.h>
#include <LcFunctionSource.h>
#include <LcGridOdePointIntegratorUpdater.h>
#include <LcLagrangeTensorElement.h>
#include <LcLinCombiner.h>
#include <LcLobattoElement1D.h>
#include <LcLorentzForceSource.h>
#include <LcLuaModuleRegistry.h>
#include <LcNodalDisContHyperUpdater.h>
#include <LcNodalDisContSrcIncrUpdater.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcPointSourceIfc.h>
#include <LcProjectOnBasisUpdater.h>
#include <LcSerendipityElement.h>
#include <LcSerendipityElement2D.h>
#include <LcSolverRegistry.h>
#include <LcWavePropagationUpdater.h>
#include <LcZeroNormalBoundaryCondition.h>
#include <LcZeroTangentBoundaryCondition.h>

#ifdef HAVE_PETSC
# include <LcFemPoissonStructUpdater.h>
#endif

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

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::GridOdePointIntegratorUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::GridOdePointIntegratorUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::GridOdePointIntegratorUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContSrcIncrUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContSrcIncrUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContSrcIncrUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::ProjectOnBasisUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::ProjectOnBasisUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::ProjectOnBasisUpdater<3> >;

#ifdef HAVE_PETSC
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FemPoissonStructUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FemPoissonStructUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::FemPoissonStructUpdater<3> >;
#endif

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EvalOnNodesUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EvalOnNodesUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::EvalOnNodesUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContHyperUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContHyperUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::NodalDisContHyperUpdater<3> >;

    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::CopyContToDisContFieldUpdater<1> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::CopyContToDisContFieldUpdater<2> >;
    new Lucee::ObjRegistry<Lucee::UpdaterIfc, Lucee::CopyContToDisContFieldUpdater<3> >;

// register boundary conditions
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::CopyBoundaryCondition>;
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::ConstBoundaryCondition>;
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::ZeroNormalBoundaryCondition>;
    new Lucee::ObjRegistry<Lucee::BoundaryCondition, Lucee::ZeroTangentBoundaryCondition>;

// register point sources
    new Lucee::ObjRegistry<Lucee::PointSourceIfc, Lucee::LorentzForceSource>;
    new Lucee::ObjRegistry<Lucee::PointSourceIfc, Lucee::CurrentSource>;
    new Lucee::ObjRegistry<Lucee::PointSourceIfc, Lucee::FunctionSource>;

// register nodal basis functions
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<1>, Lucee::LobattoElement1D>;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<2>, Lucee::SerendipityElement2D>;

    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<1>, Lucee::LagrangeTensorElement<1> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<2>, Lucee::LagrangeTensorElement<2> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<3>, Lucee::LagrangeTensorElement<3> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<4>, Lucee::LagrangeTensorElement<4> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<5>, Lucee::LagrangeTensorElement<5> >;

    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<1>, Lucee::SerendipityElement<1> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<2>, Lucee::SerendipityElement<2> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<3>, Lucee::SerendipityElement<3> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<4>, Lucee::SerendipityElement<4> >;
    new Lucee::ObjRegistry<Lucee::NodalFiniteElementIfc<5>, Lucee::SerendipityElement<5> >;


// register boundary condition library into Lucee (this needs to be
// done once here as boundary conditions are local to the slvr
// directory. Perhaps could have also done in lucee directory)
    Lucee::LuaModuleRegistry<Lucee::BoundaryCondition>::registerModule(L);

// register point source library into Lucee. Ditto comment as for
// boundary condition. See above.
    Lucee::LuaModuleRegistry<Lucee::PointSourceIfc>::registerModule(L);

// register nodal basis library into Lucee. Ditto comment as for
// boundary condition. See above.
    Lucee::LuaModuleRegistry<Lucee::NodalFiniteElementIfc<1> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::NodalFiniteElementIfc<2> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::NodalFiniteElementIfc<3> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::NodalFiniteElementIfc<4> >::registerModule(L);
    Lucee::LuaModuleRegistry<Lucee::NodalFiniteElementIfc<5> >::registerModule(L);
  }
}
