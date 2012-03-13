# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

if (PETSC_FOUND)
   include_directories(${Petsc_INCLUDE_DIRS})
endif (PETSC_FOUND)

# Create the library
set(lcslvrs_SOURCES
  LcAdvectionEquation.cpp
  LcBcUpdater.cpp
  LcBoundaryCondition.cpp
  LcConstBoundaryCondition.cpp
  LcCopyBoundaryCondition.cpp
  LcCurrentSource.cpp
  LcEdgeFaceCurlUpdater.cpp
  LcEulerEquation.cpp
  LcFaceEdgeCurlUpdater.cpp
  LcFunctionSource.cpp
  LcGridOdeIntegrator.cpp
  LcGridOdePointIntegrator.cpp
  LcGridOdePointIntegratorUpdater.cpp
  LcHyperEquation.cpp
  LcHyperEquationRegistry.cpp
  LcLinCombiner.cpp
  LcLobattoElement1D.cpp
  LcLorentzForceSource.cpp
  LcMaxwellEquation.cpp
  LcNodalFiniteElementIfc.cpp
  LcPhMaxwellEquation.cpp
  LcPointSourceIfc.cpp
  LcProjectOnBasisUpdater.cpp
  LcSerendipityElement2D.cpp
  LcSolverRegistry.cpp
  LcTenMomentEquation.cpp
  LcUpdaterIfc.cpp
  LcUpdaterStatus.cpp
  LcWavePropagationUpdater.cpp
  LcZeroNormalBoundaryCondition.cpp
  LcZeroTangentBoundaryCondition.cpp
)

if (PETSC_FOUND)
  set (lcslvrs_SOURCES ${lcslvrs_SOURCES}
    LcFemPoissonStructUpdater.cpp
  )
endif (PETSC_FOUND)

# Create the library
add_library(lcslvrs ${lcslvrs_SOURCES})
