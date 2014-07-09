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
#include <LcASquaredProjectionUpdater.h>
#include <LcATimesPUpdater.h>
#include <LcBcUpdater.h>
#include <LcBoltzmannPhiUpdater.h>
#include <LcConstBoundaryCondition.h>
#include <LcCopyBoundaryCondition.h>
#include <LcCopyContToDisContFieldUpdater.h>
#include <LcCurrentSource.h>
#include <LcEdgeFaceCurlUpdater.h>
#include <LcElectromagneticAUpdater.h>
#include <LcElectromagneticDistFuncReflectionBcUpdater.h>
#include <LcElectromagneticMomentsAtEdgesUpdater.h>
#include <LcElectromagneticZeroVelocitySource.h>
#include <LcElectromagneticZeroVelocitySourceProjection.h>
#include <LcElectrostaticPhiUpdater.h>
#include <LcEulerAxiSource.h>
#include <LcEvalOnBoundaryNodesUpdater.h>
#include <LcEvalOnNodesUpdater.h>
#include <LcFaceEdgeCurlUpdater.h>
#include <LcFieldArithmeticUpdater.h>
#include <LcFieldFunctionBoundaryCondition.h>
#include <LcFieldFunctionSource.h>
#include <LcFiveMomentNumDensityRelax.h>
#include <LcFunctionBoundaryCondition.h>
#include <LcFunctionSource.h>
#include <LcGridOdePointIntegratorUpdater.h>
#include <LcHeatFluxAtEdgeUpdater.h>
#include <LcImplicitFiveMomentSrcUpdater.h>
#include <LcImplicitTenMomentCollisionUpdater.h>
#include <LcImplicitTenMomentSrcUpdater.h>
#include <LcKineticEnergyUpdater.h>
#include <LcKineticHeatFluxAtEdgeUpdater.h>
#include <LcKineticTotalEnergyUpdater.h>
#include <LcLagrangeTensorElement.h>
#include <LcLenardBernsteinDiffUpdater.h>
#include <LcLenardBernsteinDragUpdater.h>
#include <LcLinCombiner.h>
#include <LcLobattoElement1D.h>
#include <LcLorentzForceSource.h>
#include <LcLuaModuleRegistry.h>
#include <LcMHDHamiltonianUpdater.h>
#include <LcMomentsAtEdgesUpdater.h>
#include <LcNodalCopyFaceToInteriorUpdater.h>
#include <LcNodalDgFunctionBoundaryCondition.h>
#include <LcNodalDisContHyperUpdater.h>
#include <LcNodalDisContSrcIncrUpdater.h>
#include <LcNodalFiniteElementIfc.h>
#include <LcNodalHyperDiffusionUpdater.h>
#include <LcPointSourceIfc.h>
#include <LcPositivityUpdater.h>
#include <LcPredicateUpdater.h>
#include <LcProjectOnBasisUpdater.h>
#include <LcProjectOnNodalBasisUpdater.h>
#include <LcRegisteredObjList.h>
#include <LcSOLElectronDensityInitialization.h>
#include <LcSOLElectronDensityInitialization.h>
#include <LcSOLIonDensityInitialization.h>
#include <LcSerendipityElement.h>
#include <LcSerendipityElement2D.h>
#include <LcSetPhiAtBoundaryUpdater.h>
#include <LcSolverRegistry.h>
#include <LcStairSteppedBcUpdater.h>
#include <LcTenMomLocalAnisoHeatFluxUpdater.h>
#include <LcTenMomentFluidSource.h>
#include <LcTenMomentLocalCollisionlessHeatFluxUpdater.h>
#include <LcTwoFluidMomentumRelaxSrcUpdater.h>
#include <LcVelocitiesFromMomentsUpdater.h>
#include <LcWavePropagationUpdater.h>
#include <LcZeroNormalBoundaryCondition.h>
#include <LcZeroTangentBoundaryCondition.h>

// loki includes
#include <loki/Singleton.h>

#ifdef HAVE_PETSC
# include <LcElectrostaticContPhiUpdater.h>
# include <LcFemPoissonStructUpdater.h>
# include <LcElectromagneticContAUpdater.h>
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

      .append<StairSteppedBcUpdater<1> >()
      .append<StairSteppedBcUpdater<2> >()
      .append<StairSteppedBcUpdater<3> >()
      
      // Probably used for SOL problem simulations
      .append<Lucee::BoltzmannPhiUpdater>()
      .append<Lucee::HeatFluxAtEdgeUpdater>()
      .append<Lucee::KineticEnergyUpdater>()
      .append<Lucee::KineticHeatFluxAtEdgeUpdater>()
      .append<Lucee::KineticTotalEnergyUpdater>()
      .append<Lucee::MHDHamiltonianUpdater>()
      .append<Lucee::MomentsAtEdgesUpdater>()
      .append<Lucee::ElectromagneticAUpdater>()
      .append<Lucee::ElectromagneticDistFuncReflectionBcUpdater>()
      .append<Lucee::ElectromagneticMomentsAtEdgesUpdater>()
      .append<Lucee::ElectromagneticZeroVelocitySource>()
      .append<Lucee::ElectromagneticZeroVelocitySourceProjection>()
      .append<Lucee::ElectrostaticPhiUpdater>()
      .append<Lucee::ASquaredProjectionUpdater>()
      .append<Lucee::ATimesPUpdater>()
      .append<Lucee::SetPhiAtBoundaryUpdater>()
      .append<Lucee::SOLElectronDensityInitialization>()
      .append<Lucee::SOLIonDensityInitialization>()

      // Arbitrary field arithmetic for two fields
      .append<Lucee::FieldArithmeticUpdater<1> >()
      .append<Lucee::FieldArithmeticUpdater<2> >()
      .append<Lucee::FieldArithmeticUpdater<2> >()

      .append<Lucee::GridOdePointIntegratorUpdater<1> >()
      .append<Lucee::GridOdePointIntegratorUpdater<2> >()
      .append<Lucee::GridOdePointIntegratorUpdater<3> >()

      .append<Lucee::LenardBernsteinDragUpdater>()
      .append<Lucee::LenardBernsteinDiffUpdater>()
      
      .append<Lucee::NodalCopyFaceToInteriorUpdater>()

      .append<Lucee::NodalDisContSrcIncrUpdater<1> >()
      .append<Lucee::NodalDisContSrcIncrUpdater<2> >()
      .append<Lucee::NodalDisContSrcIncrUpdater<3> >()

      .append<Lucee::NodalHyperDiffusionUpdater<1> >()
      .append<Lucee::NodalHyperDiffusionUpdater<2> >()
      .append<Lucee::NodalHyperDiffusionUpdater<3> >()

      .append<Lucee::ProjectOnBasisUpdater<1> >()
      .append<Lucee::ProjectOnBasisUpdater<2> >()
      .append<Lucee::ProjectOnBasisUpdater<3> >()

#ifdef HAVE_PETSC
      .append<Lucee::FemPoissonStructUpdater<1> >()
      .append<Lucee::FemPoissonStructUpdater<2> >()
      .append<Lucee::FemPoissonStructUpdater<3> >()

      .append<Lucee::ElectrostaticContPhiUpdater>()

      .append<Lucee::ElectromagneticContAUpdater>()
#endif

      .append<Lucee::EvalOnNodesUpdater<1> >()
      .append<Lucee::EvalOnNodesUpdater<2> >()
      .append<Lucee::EvalOnNodesUpdater<3> >()

      .append<Lucee::EvalOnBoundaryNodesUpdater<1> >()
      .append<Lucee::EvalOnBoundaryNodesUpdater<2> >()
      .append<Lucee::EvalOnBoundaryNodesUpdater<3> >()
      
      .append<Lucee::PositivityUpdater>()
    
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

      .append<Lucee::TenMomentLocalCollisionlessHeatFluxUpdater<1> >()
      .append<Lucee::TenMomentLocalCollisionlessHeatFluxUpdater<2> >()
      .append<Lucee::TenMomentLocalCollisionlessHeatFluxUpdater<3> >()

      .append<Lucee::TenMomLocalAnisoHeatFluxUpdater<1> >()
      .append<Lucee::TenMomLocalAnisoHeatFluxUpdater<2> >()
      .append<Lucee::TenMomLocalAnisoHeatFluxUpdater<3> >()

      .append<Lucee::ImplicitFiveMomentSrcUpdater<1> >()
      .append<Lucee::ImplicitFiveMomentSrcUpdater<2> >()
      .append<Lucee::ImplicitFiveMomentSrcUpdater<3> >()

      .append<Lucee::ImplicitTenMomentSrcUpdater<1> >()
      .append<Lucee::ImplicitTenMomentSrcUpdater<2> >()
      .append<Lucee::ImplicitTenMomentSrcUpdater<3> >()

      .append<TwoFluidMomentumRelaxSrcUpdater<1> >()
      .append<TwoFluidMomentumRelaxSrcUpdater<2> >()
      .append<TwoFluidMomentumRelaxSrcUpdater<3> >()

      .append<FiveMomentNumDensityRelax<1> >()
      .append<FiveMomentNumDensityRelax<2> >()
      .append<FiveMomentNumDensityRelax<3> >()

      .append<PredicateUpdater<1> >()
      .append<PredicateUpdater<2> >()
      .append<PredicateUpdater<3> >()

      .append<Lucee::VelocitiesFromMomentsUpdater>();

// register boundary conditions
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::BoundaryCondition> >
      ::Instance()
      .append<Lucee::CopyBoundaryCondition>()
      .append<Lucee::ConstBoundaryCondition>()
      .append<Lucee::ZeroNormalBoundaryCondition>()
      .append<Lucee::ZeroTangentBoundaryCondition>()
      .append<Lucee::FieldFunctionBoundaryCondition>()
      .append<Lucee::FunctionBoundaryCondition>()

      .append<Lucee::NodalDgFunctionBoundaryCondition<1> >()
      .append<Lucee::NodalDgFunctionBoundaryCondition<2> >()
      .append<Lucee::NodalDgFunctionBoundaryCondition<3> >();

// register point sources
    Loki::SingletonHolder<Lucee::RegisteredObjList<Lucee::PointSourceIfc> >
      ::Instance()
      .append<Lucee::LorentzForceSource>()
      .append<Lucee::CurrentSource>()
      .append<Lucee::FunctionSource>()
      .append<Lucee::TenMomentFluidSource>()
      .append<Lucee::FieldFunctionSource>()
      .append<Lucee::EulerAxiSource>();

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
