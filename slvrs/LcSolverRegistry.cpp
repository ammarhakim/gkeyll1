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
#include <LcFieldFunctionBoundaryCondition.h>
#include <LcFunctionBoundaryCondition.h>
#include <LcFunctionSource.h>
#include <LcGridOdePointIntegratorUpdater.h>
#include <LcImplicitFiveMomentSrcUpdater.h>
#include <LcImplicitTenMomentCollisionUpdater.h>
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
#include <LcProjectOnNodalBasisUpdater.h>
#include <LcRegisteredObjList.h>
#include <LcSerendipityElement.h>
#include <LcSerendipityElement2D.h>
#include <LcSolverRegistry.h>
#include <LcTenMomentFluidSource.h>
#include <LcWavePropagationUpdater.h>
#include <LcZeroNormalBoundaryCondition.h>
#include <LcZeroTangentBoundaryCondition.h>

// loki includes
#include <loki/Singleton.h>

#ifdef HAVE_PETSC
# include <LcFemPoissonStructUpdater.h>
#endif

namespace Lucee
{
  void
  registerSolverObjects(Lucee::LuaState& L)
  {
// register updaters
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::UpdaterIfc> >
      ::Instance()
      .append<Lucee::LinCombiner<1> >()
      .append<Lucee::LinCombiner<2> >()
      .append<Lucee::LinCombiner<3> >()

      .append<Lucee::WavePropagationUpdater<1> >()
      .append<Lucee::WavePropagationUpdater<2> >()
      .append<Lucee::WavePropagationUpdater<3> >()

      .append<Lucee::FaceEdgeCurlUpdater<1> >()
      .append<Lucee::FaceEdgeCurlUpdater<2> >()
      .append<Lucee::FaceEdgeCurlUpdater<3> >()

      .append<Lucee::EdgeFaceCurlUpdater<1> >()
      .append<Lucee::EdgeFaceCurlUpdater<2> >()
      .append<Lucee::EdgeFaceCurlUpdater<3> >()

      .append<Lucee::BcUpdater<1> >()
      .append<Lucee::BcUpdater<2> >()
      .append<Lucee::BcUpdater<3> >()

      .append<Lucee::GridOdePointIntegratorUpdater<1> >()
      .append<Lucee::GridOdePointIntegratorUpdater<2> >()
      .append<Lucee::GridOdePointIntegratorUpdater<3> >()

      .append<Lucee::NodalDisContSrcIncrUpdater<1> >()
      .append<Lucee::NodalDisContSrcIncrUpdater<2> >()
      .append<Lucee::NodalDisContSrcIncrUpdater<3> >()

      .append<Lucee::ProjectOnBasisUpdater<1> >()
      .append<Lucee::ProjectOnBasisUpdater<2> >()
      .append<Lucee::ProjectOnBasisUpdater<3> >()

#ifdef HAVE_PETSC
      .append<Lucee::FemPoissonStructUpdater<1> >()
      .append<Lucee::FemPoissonStructUpdater<2> >()
      .append<Lucee::FemPoissonStructUpdater<3> >()
#endif

      .append<Lucee::EvalOnNodesUpdater<1> >()
      .append<Lucee::EvalOnNodesUpdater<2> >()
      .append<Lucee::EvalOnNodesUpdater<3> >()
    
      .append<Lucee::ProjectOnNodalBasisUpdater<1> >()
      .append<Lucee::ProjectOnNodalBasisUpdater<2> >()
      .append<Lucee::ProjectOnNodalBasisUpdater<3> >()

      .append<Lucee::NodalDisContHyperUpdater<1> >()
      .append<Lucee::NodalDisContHyperUpdater<2> >()
      .append<Lucee::NodalDisContHyperUpdater<3> >()

      .append<Lucee::CopyContToDisContFieldUpdater<1> >()
      .append<Lucee::CopyContToDisContFieldUpdater<2> >()
      .append<Lucee::CopyContToDisContFieldUpdater<3> >()

      .append<Lucee::ImplicitTenMomentCollisionUpdater<1> >()
      .append<Lucee::ImplicitTenMomentCollisionUpdater<2> >()
      .append<Lucee::ImplicitTenMomentCollisionUpdater<3> >()

      .append<Lucee::ImplicitFiveMomentSrcUpdater<1> >()
      .append<Lucee::ImplicitFiveMomentSrcUpdater<2> >()
      .append<Lucee::ImplicitFiveMomentSrcUpdater<3> >();

// register boundary conditions
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::BoundaryCondition> >
      ::Instance()
      .append<Lucee::CopyBoundaryCondition>()
      .append<Lucee::ConstBoundaryCondition>()
      .append<Lucee::ZeroNormalBoundaryCondition>()
      .append<Lucee::ZeroTangentBoundaryCondition>()
      .append<Lucee::FieldFunctionBoundaryCondition>()
      .append<Lucee::FunctionBoundaryCondition>();

// register point sources
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::PointSourceIfc> >
      ::Instance()
      .append<Lucee::LorentzForceSource>()
      .append<Lucee::CurrentSource>()
      .append<Lucee::FunctionSource>()
      .append<Lucee::TenMomentFluidSource>();

// register nodal basis functions
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::NodalFiniteElementIfc<1> > >
      ::Instance()
      .append<Lucee::LobattoElement1D>()
      .append<Lucee::LagrangeTensorElement<1> >()
      .append<Lucee::SerendipityElement<1> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::NodalFiniteElementIfc<2> > >
      ::Instance()
      .append<Lucee::SerendipityElement2D>()
      .append<Lucee::LagrangeTensorElement<2> >()
      .append<Lucee::SerendipityElement<2> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::NodalFiniteElementIfc<3> > >
      ::Instance()
      .append<Lucee::LagrangeTensorElement<3> >()
      .append<Lucee::SerendipityElement<3> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::NodalFiniteElementIfc<4> > >
      ::Instance()
      .append<Lucee::LagrangeTensorElement<4> >()
      .append< Lucee::SerendipityElement<4> >();

    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::NodalFiniteElementIfc<5> > >
      ::Instance()
      .append<Lucee::LagrangeTensorElement<5> >()
      .append<Lucee::SerendipityElement<5> >();

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
