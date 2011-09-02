# Project information
project(lucee)

set(VERSION_MAJOR "0")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")
set(LUCEE_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8)

# initialize build
INCLUDE(${CMAKE_SOURCE_DIR}/CMake/TxInit.cmake)
if (WIN32)
  # do nothing to find Fortran
else (WIN32)
  # determine Fortran compiler
  include(CMakeDetermineFortranCompiler)
  # now initialize Fortran stuff
  include (${CMAKE_SOURCE_DIR}/CMake/TxFortranChecks.cmake)
endif (WIN32)

######################################################################
#
# Set for inclusion in config.h
#
######################################################################

set(LUCEE_REV ${PROJECT_REV})
set(CXX ${CMAKE_CXX_COMPILER})
# MESSAGE("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
set(BUILD_FLAGS_VAR  CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
# MESSAGE("BUILD_FLAGS_VAR = ${BUILD_FLAGS_VAR}")
set(BUILD_FLAGS_VAL "${${BUILD_FLAGS_VAR}}")
# MESSAGE("BUILD_FLAGS_VAL = ${BUILD_FLAGS_VAL}")
set(CXXFLAGS "${BUILD_FLAGS_VAL} ${CMAKE_CXX_FLAGS}")
# MESSAGE("CXXFLAGS = ${CXXFLAGS}")
# SET(HAVE_CMAKE CACHE STRING "Whether built with CMake")
set(HAVE_CMAKE 1 CACHE STRING "Whether built with CMake")

######################################################################
#
# Find packages
#
######################################################################

# Prerequisites for finding packages
include(${CMAKE_SOURCE_DIR}/CMake/Modules/TxFindPackage.cmake)

# find txbase
find_package(TxBase REQUIRED)
if (TXBASE_FOUND)
   include_directories(${TxBase_INCLUDE_DIRS})
endif (TXBASE_FOUND)

# find HDF5
find_package(Hdf5 REQUIRED)
if (HDF5_FOUND)
   include_directories(${Hdf5_INCLUDE_DIRS})
endif (HDF5_FOUND)
if (NOT WIN32)
  find_package(Z REQUIRED)
endif ()

# find GSL
find_package(Gsl REQUIRED)
if (GSL_FOUND)
   include_directories(${Gsl_INCLUDE_DIRS})
endif (GSL_FOUND)

# find FFTW3
find_package(Fftw3)

IF(APPLE)
  # add the Accelerate framework on Mac OS X
  FIND_LIBRARY(ACCELERATE_LIB Accelerate)
  MARK_AS_ADVANCED(ACCELERATE_LIB)
  SET(EXTRA_LIBS ${ACCELERATE_LIB})
elseif(WIN32)
  # CLapack is used on Windows
  find_package(CLapackCMake REQUIRED)
  set(EXTRA_LIBS ${CLapackCMake_LIBRARIES})
else(APPLE)
  # find LAPACK
  find_package(LAPACK REQUIRED)
  if (LAPACK_FOUND)
    add_definitions(-DHAVE_LAPACK)
    SET(EXTRA_LIBS ${LAPACK_LIBRARIES})
  endif (LAPACK_FOUND)
endif(APPLE)

if (ENABLE_PARALLEL AND MPI_INCLUDE_PATH AND NOT TX_HAVE_MPI_COMPILER_WRAPPER)
  include_directories("${MPI_INCLUDE_PATH}")
endif ()

# Create configure file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in
   ${CMAKE_CURRENT_BINARY_DIR}/config.h)

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

# lcgeomprim
add_test(lcgeomprim unit/lcgeomprim)

# lcstructgridfield
add_test(lcstructgridfield unit/lcstructgridfield)
