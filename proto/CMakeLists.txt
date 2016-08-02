# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/proto)
include_directories(${PROJECT_SOURCE_DIR}/etc/quadrule)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)

if (FFTW3_FOUND)
   include_directories(${Fftw3_INCLUDE_DIRS})
endif (FFTW3_FOUND)

if (PETSC_FOUND)
   include_directories(${Petsc_INCLUDE_DIRS})
endif (PETSC_FOUND)

set (proto_SOURCES
  LcBGKCollUpdater.cpp	
  LcBoltzmannCollisionOperator.cpp
  LcCenterOfMassUpdater.cpp
  LcConstGravitySrcUpdater.cpp
  LcConstructLinearOperatorMatrix.cpp
  LcCopy1DTo2DNodalField.cpp
  LcCopyNodalFields.cpp
  LcDGDiffusionUpdater1D.cpp
  LcDistFuncMomentCalc1D.cpp
  LcDistFuncMomentCalc1DFrom3D.cpp
  LcDistFuncMomentCalc1DFrom4D.cpp
  LcDistFuncMomentCalc2D.cpp
  LcDistFuncMomentCalc3D.cpp
  LcDistFuncMomentCalcCDIMFromVDIM.cpp
  LcDistFuncMomentCalcWeighted2D.cpp
  LcDistFuncMomentCalcWeighted3D.cpp
  LcDistFuncReflectionBcUpdater.cpp
  LcDualSymmetry1D2VUpdater.cpp
  LcETGAdiabaticPotentialUpdater.cpp
  LcETGAdiabaticPotentialUpdater3D.cpp
  LcETGAdjointSource.cpp
  LcETGFreeEnergy.cpp
  LcETGInitializeDensity.cpp
  LcETGInitializeDensity5D.cpp
  LcETGInnerProduct.cpp
  LcEnergyFromStreamAndVortUpdater.cpp
  LcEnergyFromStreamFunctionUpdater.cpp
  LcEnstrophyUpdater.cpp
  LcFiniteVolumeToLinearDGUpdater.cpp
  LcGradLinearDGToFiniteVolumeUpdater.cpp
  LcInitNodesFromMatrixMarketUpdater.cpp
  LcIntegrateField.cpp
  LcIntegrateFieldAlongLine.cpp
  LcIntegrateFieldProduct.cpp
  LcIntegrateGeneralField.cpp
  LcIntegrateNodalField.cpp
  LcIonizationSource.cpp
  LcLinEmGke1dHamilPertUpdater.cpp
  LcMaxwellDistInit.cpp
  LcMaxwellTm2DUpdater.cpp
  LcModalDg1DDiffusionUpdater.cpp
  LcModalDg1DHyperDiffusionUpdater.cpp
  LcModalDg1DLocalDGUpdater.cpp
  LcModalDg1DSymmetricDDGUpdater.cpp
  LcModalDg1DUpdater.cpp
  LcModalDgLimiter1DUpdater.cpp
  LcModalL2NormUpdater.cpp
  LcMusclHancock1DUpdater.cpp
  LcNeutralDragForceSource.cpp
  LcNodalCopy1DTo3DFieldUpdater.cpp
  LcNodalCopy2DTo4DFieldUpdater.cpp
  LcNodalCopy3DTo5DFieldUpdater.cpp
  LcNodalDgConstGravitySrcUpdater.cpp
  LcNodalGradientUpdater.cpp
  LcNodalPoissonBracketUpdater.cpp
  LcNonLinEmGke1dHamilUpdater.cpp
  LcNormGradPhiUpdater.cpp
  LcPoissonBracketAdvectionEquation1D.cpp
  LcPoissonBracketAdvectionEquation4D.cpp
  LcPoissonBracketAdvectionEquation5D.cpp
  LcPoissonBracketCanonical2D.cpp
  LcPoissonBracketCanonical3D.cpp
  LcPoissonBracketCanonical4D.cpp
  LcPoissonBracketEquation.cpp
  LcPoissonBracketEquationRegistry.cpp
  LcPoissonBracketGyroEquation4D.cpp
  LcPoissonBracketGyroEquation5D.cpp
  LcPoissonBracketImpUpdater.cpp
  LcPoissonBracketOptUpdater.cpp
  LcPoissonBracketSOL3D.cpp
  LcPoissonBracketSimpleUpdater.cpp
  LcPoissonBracketUpdater.cpp
  LcPoissonBracketVlasovPoissonEquation3D.cpp
  LcProtoSolverRegistry.cpp
  LcRecordFieldDerivInCell.cpp
  LcRecordFieldInCell.cpp
  LcRecoveryDG1DUpdater.cpp
  LcRecoveryDG3DUpdater.cpp
  LcRectFourthOrderCentralDiffUpdater.cpp
  LcRectSecondOrderCentralDiffUpdater.cpp
  LcRunningAverageOfFieldCalc.cpp
  LcSOLFluxAcrossEdgeCalc.cpp
  LcSOLSetPotentialAtBoundary.cpp
  LcSOLTotalEnergyCalc.cpp
  LcSOLTotalIntegralCalc.cpp
  LcSOLUpperXPotentialBcUpdater.cpp
  LcSetSingleNodeToOneUpdater.cpp
  LcSheathParticleSource1x1v.cpp
  LcSimpleSmoothToC0Updater.cpp
  LcSmoothQuadPhiToC1Updater.cpp
  LcThreeWaveInteractModSrcUpdater.cpp
  LcThreeWaveInteractSrcUpdater.cpp
  LcZonalAverageCalc3D.cpp
  LcZonalVelocity1DUpdater.cpp
)

if (FFTW3_FOUND)
  set (proto_SOURCES ${proto_SOURCES} LcPeriodicPoisson2DUpdater.cpp
  LcPeriodicCollisionlessHeatFluxUpdater.cpp )
  if (ENABLE_PARALLEL)
    set (proto_SOURCES ${proto_SOURCES} LcPeriodicParallelPoisson2DUpdater.cpp)
  endif (ENABLE_PARALLEL)
endif (FFTW3_FOUND)

if (PETSC_FOUND)
  set (proto_SOURCES ${proto_SOURCES}
    LcContFromDisContUpdater.cpp
  )
endif (PETSC_FOUND)

# Create the library
add_library(lcproto ${proto_SOURCES})

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY ./
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
  )
  
  install(TARGETS lcproto
    ARCHIVE DESTINATION lib
  )
endif ()
