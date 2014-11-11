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
  LcCenterOfMassUpdater.cpp
  LcConstGravitySrcUpdater.cpp
  LcCopy1DTo2DNodalField.cpp
  LcDGDiffusionUpdater1D.cpp
  LcDistFuncMomentCalc1D.cpp
  LcDistFuncMomentCalc2D.cpp
  LcDistFuncMomentCalcWeighted2D.cpp
  LcDistFuncReflectionBcUpdater.cpp
  LcETGAdiabaticPotentialUpdater.cpp
  LcETGInitializeDensity.cpp
  LcEnergyFromStreamAndVortUpdater.cpp
  LcEnergyFromStreamFunctionUpdater.cpp
  LcEnstrophyUpdater.cpp
  LcIntegrateField.cpp
  LcIntegrateFieldAlongLine.cpp
  LcIntegrateGeneralField.cpp
  LcIntegrateNodalField.cpp
  LcLinEmGke1dHamilPertUpdater.cpp
  LcMaxwellTm2DUpdater.cpp
  LcModalDg1DDiffusionUpdater.cpp
  LcModalDg1DHyperDiffusionUpdater.cpp
  LcModalDg1DLocalDGUpdater.cpp
  LcModalDg1DSymmetricDDGUpdater.cpp
  LcModalDg1DUpdater.cpp
  LcModalDgLimiter1DUpdater.cpp
  LcModalL2NormUpdater.cpp
  LcMusclHancock1DUpdater.cpp
  LcNodalCopy2DTo4DFieldUpdater.cpp
  LcNodalDgConstGravitySrcUpdater.cpp
  LcNodalGradientUpdater.cpp
  LcNodalPoissonBracketUpdater.cpp
  LcNonLinEmGke1dHamilUpdater.cpp
  LcNormGradPhiUpdater.cpp
  LcPoissonBracketAdvectionEquation4D.cpp
  LcPoissonBracketCanonical2D.cpp
  LcPoissonBracketCanonical4D.cpp
  LcPoissonBracketEquation.cpp
  LcPoissonBracketEquationRegistry.cpp
  LcPoissonBracketGyroEquation4D.cpp
  LcPoissonBracketUpdater.cpp
  LcProtoSolverRegistry.cpp
  LcRecordFieldDerivInCell.cpp
  LcRecordFieldInCell.cpp
  LcRectSecondOrderCentralDiffUpdater.cpp
  LcSheathParticleSource1x1v.cpp
  LcSimpleSmoothToC0Updater.cpp
  LcSmoothQuadPhiToC1Updater.cpp
  LcThreeWaveInteractModSrcUpdater.cpp
  LcThreeWaveInteractSrcUpdater.cpp
)

if (FFTW3_FOUND)
  set (proto_SOURCES ${proto_SOURCES} LcPeriodicPoisson2DUpdater.cpp)
endif (FFTW3_FOUND)

if (PETSC_FOUND)
  set (proto_SOURCES ${proto_SOURCES}
    LcContFromDisContUpdater.cpp
  )
endif (PETSC_FOUND)

# Create the library
add_library(lcproto ${proto_SOURCES})


