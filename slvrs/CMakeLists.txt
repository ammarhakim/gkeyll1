# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/etc/quadrule)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

if (PETSC_FOUND)
   include_directories(${Petsc_INCLUDE_DIRS})
endif (PETSC_FOUND)

# Create the library
set(lcslvrs_SOURCES
  LcASquaredProjectionUpdater.cpp
  LcATimesPUpdater.cpp
  LcAdvectionEquation.cpp
  LcAdvectionEquationTemplated.cpp
  LcAuxAdvectionEquation.cpp
  LcBcUpdater.cpp
  LcBoltzmannPhiUpdater.cpp
  LcBoundaryCondition.cpp
  LcCompletePolynomialElement.cpp
  LcConstBoundaryCondition.cpp
  LcCopyBoundaryCondition.cpp
  LcCopyContToDisContFieldUpdater.cpp
  LcCurrentSource.cpp
  LcDivEquation.cpp
  LcEdgeFaceCurlUpdater.cpp
  LcElectromagneticAUpdater.cpp
  LcElectromagneticDistFuncReflectionBcUpdater.cpp
  LcElectromagneticMomentsAtEdgesUpdater.cpp
  LcElectromagneticZeroVelocitySource.cpp
  LcElectromagneticZeroVelocitySourceProjection.cpp
  LcElectrostaticPhiUpdater.cpp
  LcEulerAxiSource.cpp
  LcEulerEquation.cpp
  LcEvalOnBoundaryNodesUpdater.cpp
  LcEvalOnCentroidUpdater.cpp
  LcEvalOnNodesUpdater.cpp
  LcFaceEdgeCurlUpdater.cpp
  LcFieldArithmeticUpdater.cpp
  LcFieldFunctionBoundaryCondition.cpp
  LcFieldFunctionSource.cpp
  LcFiveMomentNumDensityRelax.cpp
  LcFunctionBoundaryCondition.cpp
  LcFunctionSource.cpp
  LcGradEquation.cpp
  LcGridOdeIntegrator.cpp
  LcGridOdePointIntegrator.cpp
  LcGridOdePointIntegratorUpdater.cpp
  LcHeatFluxAtEdge3DUpdater.cpp
  LcHeatFluxAtEdgeUpdater.cpp
  LcHyperEquation.cpp
  LcHyperEquationRegistry.cpp
  LcImplicitFiveMomentSrcUpdater.cpp
  LcImplicitTenMomentCollisionUpdater.cpp
  LcImplicitTenMomentSrcUpdater.cpp
  LcImplicitTwentyMomentSrcUpdater.cpp
  LcImplicitTwentyMomentCollisionUpdater.cpp
  LcIsoThermEulerEquation.cpp
  LcKineticEnergyUpdater.cpp
  LcKineticHeatFluxAtEdge3DUpdater.cpp
  LcKineticHeatFluxAtEdgeUpdater.cpp
  LcKineticTotalEnergyUpdater.cpp
  LcLagrangeTensorBasisCalc.cpp
  LcLagrangeTensorElement.cpp
  LcLenardBernsteinDiff3DUpdater.cpp
  LcLenardBernsteinDiffUpdater.cpp
  LcLenardBernsteinDrag3DUpdater.cpp
  LcLenardBernsteinDragUpdater.cpp
  LcLenardBernsteinVParEquation.cpp
  LcLinCombiner.cpp
  LcLobattoElement1D.cpp
  LcLorentzForceSource.cpp
  LcMHDHamiltonianUpdater.cpp
  LcMaxwellEquation.cpp
  LcMomentsAtEdges3DUpdater.cpp
  LcMomentsAtEdgesUpdater.cpp
  LcNodalCopyFaceToInteriorUpdater.cpp
  LcNodalDgCopyBoundaryCondition.cpp
  LcNodalDgFunctionBoundaryCondition.cpp
  LcNodalDgZeroNormalBoundaryCondition.cpp
  LcNodalDgZeroTangentBoundaryCondition.cpp
  LcNodalDisContHyperUpdater.cpp
  LcNodalDisContSrcIncrUpdater.cpp
  LcNodalFiniteElementIfc.cpp
  LcNodalHyperDiffusionUpdater.cpp
  LcNodalPositiveFilterUpdater.cpp
  LcNodalVlasovUpdater.cpp
  LcPhMaxwellEquation.cpp
  LcPointSourceIfc.cpp
  LcPositivityUpdater.cpp
  LcPredicateUpdater.cpp
  LcProjectOnBasisUpdater.cpp
  LcProjectOnNodalBasisUpdater.cpp
  LcReflectionBoundaryCondition.cpp
  LcSOL3DElectronTempAtWallCalc.cpp
  LcSOL3DElectrostaticDistFuncReflectionBCUpdater.cpp
  LcSOLElectronDensityInitialization.cpp
  LcSOLElectronDensityInitialization.cpp
  LcSOLIonDensityInitialization.cpp
  LcSerendipityElement.cpp
  LcSerendipityElement2D.cpp
  LcSerendipityElement2DDiffusionOutput
  LcSerendipityElementDiffusionOutput
  LcSerendipityElementHyperDiffusionOutput
  LcSetPhiAtBoundaryUpdater.cpp
  LcSolverRegistry.cpp
  LcStairSteppedBcUpdater.cpp
  LcRadialStairSteppedBcUpdater.cpp
  LcTenMomLocalAnisoHeatFluxUpdater.cpp
  LcTenMomentEquation.cpp
  LcTenMomentFluidSource.cpp
  LcTwentyMomentEquation.cpp
  LcTwentyMomentFluidSource.cpp
  LcTenMomentLocalCollisionlessHeatFluxUpdater.cpp
  LcTwoFluidMomentumRelaxSrcUpdater.cpp
  LcUpdaterIfc.cpp
  LcUpdaterStatus.cpp
  LcVelocitiesFromMoments3DUpdater.cpp
  LcVelocitiesFromMomentsUpdater.cpp
  LcWavePropagationUpdater.cpp
  LcZeroNormalBoundaryCondition.cpp
  LcZeroTangentBoundaryCondition.cpp
)

if (PETSC_FOUND)
  set (lcslvrs_SOURCES ${lcslvrs_SOURCES}
    LcElectromagneticContAUpdater.cpp
    LcElectrostaticContPhiUpdater.cpp
    LcFemGKPoissonStructUpdater.cpp 
    LcFemPoissonStructUpdater.cpp
  )
endif (PETSC_FOUND)

# Create the library
add_library(lcslvrs ${lcslvrs_SOURCES})

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY ./
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
  )
  
  install(TARGETS lcslvrs
    ARCHIVE DESTINATION lib
  )
endif ()
