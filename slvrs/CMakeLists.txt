# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)


# Create the library
add_library(lcslvrs
  LcAdvectionEquation.cpp
  LcBcUpdater.cpp
  LcBoundaryCondition.cpp
  LcCopyBoundaryCondition.cpp
  LcEdgeFaceCurlUpdater.cpp
  LcEulerEquation.cpp
  LcFaceEdgeCurlUpdater.cpp
  LcGridOdeIntegrator.cpp
  LcGridOdePointIntegrator.cpp
  LcHyperEquation.cpp
  LcHyperEquationRegistry.cpp
  LcLinCombiner.cpp
  LcLorentzForceSource.cpp
  LcMaxwellEquation.cpp
  LcPhMaxwellEquation.cpp
  LcPointSourceIfc.cpp
  LcSolverRegistry.cpp
  LcUpdaterIfc.cpp
  LcUpdaterStatus.cpp
  LcWavePropagationUpdater.cpp
  LcZeroNormalBoundaryCondition.cpp
  LcZeroTangentBoundaryCondition.cpp
)
