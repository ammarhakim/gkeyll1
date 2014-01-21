# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/pic)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)

if (FFTW3_FOUND)
   include_directories(${Fftw3_INCLUDE_DIRS})
endif (FFTW3_FOUND)

if (PETSC_FOUND)
   include_directories(${Petsc_INCLUDE_DIRS})
endif (PETSC_FOUND)

# Create the library
add_library(lcpic
  LcParticleSpecies.cpp
)


