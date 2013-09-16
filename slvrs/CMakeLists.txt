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
  LcAdvectionEquation.cpp
  LcAuxAdvectionEquation.cpp
  LcBcUpdater.cpp
  LcBoltzmannPhiUpdater.cpp
  LcBoundaryCondition.cpp
  LcConstBoundaryCondition.cpp
  LcCopyBoundaryCondition.cpp
  LcCopyContToDisContFieldUpdater.cpp
  LcCurrentSource.cpp
  LcDivEquation.cpp
  LcEdgeFaceCurlUpdater.cpp
  LcEulerEquation.cpp
  LcEvalOnBoundaryNodesUpdater.cpp
  LcEvalOnNodesUpdater.cpp
  LcFaceEdgeCurlUpdater.cpp
  LcFieldFunctionBoundaryCondition.cpp
  LcFunctionBoundaryCondition.cpp
  LcFunctionSource.cpp
  LcGradEquation.cpp
  LcGridOdeIntegrator.cpp
  LcGridOdePointIntegrator.cpp
  LcGridOdePointIntegratorUpdater.cpp
  LcHeatFluxAtEdgeUpdater.cpp
  LcHyperEquation.cpp
  LcHyperEquationRegistry.cpp
  LcImplicitFiveMomentSrcUpdater.cpp
  LcImplicitTenMomentCollisionUpdater.cpp
  LcImplicitTenMomentSrcUpdater.cpp
  LcKineticHeatFluxAtEdgeUpdater.cpp
  LcKineticTotalEnergyUpdater.cpp
  LcLagrangeTensorBasisCalc.cpp
  LcLagrangeTensorElement.cpp
  LcLenardBernsteinDiffUpdater.cpp
  LcLenardBernsteinDragUpdater.cpp
  LcLenardBernsteinVParEquation.cpp
  LcLinCombiner.cpp
  LcLobattoElement1D.cpp
  LcLorentzForceSource.cpp
  LcMaxwellEquation.cpp
  LcNodalCopyFaceToInteriorUpdater.cpp
  LcNodalDgFunctionBoundaryCondition.cpp
  LcNodalDisContHyperUpdater.cpp
  LcNodalDisContSrcIncrUpdater.cpp
  LcNodalFiniteElementIfc.cpp
  LcNodalHyperDiffusionUpdater.cpp
  LcPhMaxwellEquation.cpp
  LcPointSourceIfc.cpp
  LcPositivityUpdater.cpp
  LcProjectOnBasisUpdater.cpp
  LcProjectOnNodalBasisUpdater.cpp
  LcSerendipityElement.cpp
  LcSerendipityElement2D.cpp
  LcSerendipityElement2DDiffusionOutput
  LcSerendipityElementDiffusionOutput
  LcSerendipityElementHyperDiffusionOutput
  LcSolverRegistry.cpp
  LcTenMomentEquation.cpp
  LcTenMomentFluidSource.cpp
  LcTenMomentLocalCollisionlessHeatFluxUpdater.cpp
  LcUpdaterIfc.cpp
  LcUpdaterStatus.cpp
  LcVelocitiesFromMomentsUpdater.cpp
  LcWavePropagationUpdater.cpp
  LcZeroNormalBoundaryCondition.cpp
  LcZeroTangentBoundaryCondition.cpp
)

if (PETSC_FOUND)
  set (lcslvrs_SOURCES ${lcslvrs_SOURCES}
    LcFemPoissonStructUpdater.cpp
    LcElectrostaticContPhiUpdater.cpp
    LcElectrostaticPhiUpdater.cpp
  )
endif (PETSC_FOUND)

# Create the library
add_library(lcslvrs ${lcslvrs_SOURCES})
