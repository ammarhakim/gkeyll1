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
  LcCopy1DTo2DNodalField.cpp
  LcDGDiffusionUpdater1D.cpp
  LcDistFuncMomentCalc1D.cpp
  LcDistFuncReflectionBcUpdater.cpp
  LcEnergyFromStreamAndVortUpdater.cpp
  LcEnergyFromStreamFunctionUpdater.cpp
  LcEnstrophyUpdater.cpp
  LcIntegrateFieldAlongLine.cpp
  LcIntegrateNodalField.cpp
  LcMaxwellTm2DUpdater.cpp
  LcModalDg1DDiffusionUpdater.cpp
  LcModalDg1DHyperDiffusionUpdater.cpp
  LcModalDg1DLocalDGUpdater.cpp
  LcModalDg1DSymmetricDDGUpdater.cpp
  LcModalDg1DUpdater.cpp
  LcModalDgLimiter1DUpdater.cpp
  LcModalL2NormUpdater.cpp
  LcMusclHancock1DUpdater.cpp
  LcNodalGradientUpdater.cpp
  LcNodalPoissonBracketUpdater.cpp
  LcNormGradPhiUpdater.cpp
  LcProtoSolverRegistry.cpp
  LcRecordFieldDerivInCell.cpp
  LcRecordFieldInCell.cpp
  LcRectSecondOrderCentralDiffUpdater.cpp
  LcSheathParticleSource1x1v.cpp
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


