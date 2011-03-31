# Project information
project(lucee)

set(VERSION_MAJOR "0")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8)

# initialize build
INCLUDE(${CMAKE_SOURCE_DIR}/CMake/TxInit.cmake)

# determine Fortran compiler
include(CMakeDetermineFortranCompiler)
# now initialize Fortran stuff
include (${CMAKE_SOURCE_DIR}/CMake/TxFortranChecks.cmake)

# find packages
# Prerequisites for finding packages
include(${CMAKE_SOURCE_DIR}/CMake/Modules/TxFindPackage.cmake)

# find HDF5
find_package(Hdf5 REQUIRED)
if (HDF5_FOUND)
   add_definitions(-DHAVE_HDF5)
   add_definitions(-DNEW_H5S_SELECT_HYPERSLAB_IFC)
   include_directories(${Hdf5_INCLUDE_DIRS})
endif (HDF5_FOUND)

# find GSL
find_package(Gsl REQUIRED)
if (GSL_FOUND)
   add_definitions(-DHAVE_GSL)
   include_directories(${Gsl_INCLUDE_DIRS})
endif (GSL_FOUND)

IF(APPLE)
  # add the Accelerate framework on Mac OS X
  FIND_LIBRARY(ACCELERATE_LIB Accelerate)
  MARK_AS_ADVANCED(ACCELERATE_LIB)
  SET(EXTRA_LIBS ${ACCELERATE_LIB})
else(APPLE)
  # find LAPACK
  find_package(LAPACK REQUIRED)
  if (LAPACK_FOUND)
    add_definitions(-DHAVE_LAPACK)
    SET(EXTRA_LIBS ${LAPACK_LIBRARIES})
  endif (LAPACK_FOUND)
endif(APPLE)

# Create configure file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

install(FILES svninfo.txt DESTINATION include)

# Build code in loki
add_subdirectory(etc/loki/src)

# Build code in lua
add_subdirectory(etc/lua/src)

# Build code in lib
add_subdirectory(lib)

# Build code in grids
add_subdirectory(grids)

# Build code in datastructs
add_subdirectory(datastructs)

# Build code in slvrs
add_subdirectory(slvrs)

# Build code in proto
add_subdirectory(proto)

# Build library in rte directory
add_subdirectory(rte)

# Build Lucee executable
add_subdirectory(lucee)

# Build all unit tests
add_subdirectory(unit)

# Add individual tests
ENABLE_TESTING()

# lcarray
add_test(lcarray unit/lcarray)

# lccolmajorindexer
add_test(lccolmajorindexer unit/lccolmajorindexer)

# lcrowmajorindexer
add_test(lcrowmajorindexer unit/lcrowmajorindexer)

# lcfixedvector
add_test(lcfixedvector unit/lcfixedvector)

# lcvector
add_test(lcvector unit/lcvector)

# lcmatrix
add_test(lcmatrix unit/lcmatrix)

# lccmdlineargs
add_test(lccmdlineargs unit/lccmdlineargs)

# lcregion
add_test(lcregion unit/lcregion)

# lcfield
add_test(lcfield unit/lcfield)

# lcrowmajorsequencer
add_test(lcrowmajorsequencer unit/lcrowmajorsequencer)

# lccolmajorsequencer
add_test(lccolmajorsequencer unit/lccolmajorsequencer)

# lcdirsequencer
add_test(lcdirsequencer unit/lcdirsequencer)

# lcluatable
add_test(lcluatable unit/lcluatable)

# lcluamodule
add_test(lcluamodule unit/lcluamodule)

# lclogger
add_test(lclogger unit/lclogger)

# lchdf5io
add_test(lchdf5io unit/lchdf5io)

# lcvec3
add_test(lcvec3 unit/lcvec3)

# lcrectcartgrid
add_test(lcrectcartgrid unit/lcrectcartgrid)

# lcluafunction
add_test(lcluafunction unit/lcluafunction)

# lclincombiner
add_test(lclincombiner unit/lclincombiner)

# lcfieldptr
add_test(lcfieldptr unit/lcfieldptr)

# lcrectcoordsys
add_test(lcrectcoordsys unit/lcrectcoordsys)

# lctabledescription
add_test(lctabledescription unit/lctabledescription)

# lcluahelp
add_test(lcluahelp unit/lcluahelp)
