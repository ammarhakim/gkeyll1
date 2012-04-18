# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/proto)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

if (FFTW3_FOUND)
   include_directories(${Fftw3_INCLUDE_DIRS})
endif (FFTW3_FOUND)

if (PETSC_FOUND)
   include_directories(${Petsc_INCLUDE_DIRS})
endif (PETSC_FOUND)

set (proto_SOURCES
  LcEnergyFromStreamAndVortUpdater.cpp
  LcEnergyFromStreamFunctionUpdater.cpp
  LcEnstrophyUpdater.cpp
  LcMaxwellTm2DUpdater.cpp
  LcModalDg1DUpdater.cpp
  LcModalDgLimiter1DUpdater.cpp
  LcMusclHancock1DUpdater.cpp
  LcNodalGradientUpdater.cpp
  LcNodalPoissonBracketUpdater.cpp
  LcProtoSolverRegistry.cpp
  LcRectSecondOrderCentralDiffUpdater.cpp
)

if (FFTW3_FOUND)
  set (proto_SOURCES ${proto_SOURCES} LcPeriodicPoisson2DUpdater.cpp)
endif (FFTW3_FOUND)

# Create the library
add_library(lcproto ${proto_SOURCES})


