# Project information
project(lucee)

set(VERSION_MAJOR "0")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")
set(LUCEE_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# Required version
cmake_minimum_required(VERSION 2.8.6)

######################################################################
#
# Usual initialization stuff
#
######################################################################

if (WIN32)
  set(NOFORTRAN TRUE)
endif ()
include(${CMAKE_SOURCE_DIR}/CMake/TxInit.cmake)
set(NAUTILUS_REV "${PROJECT_REV}")

######################################################################
#
# Set for inclusion in config.h
#
######################################################################

set(CXX ${CMAKE_CXX_COMPILER})
set(BUILD_FLAGS_VAR  CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
set(BUILD_FLAGS_VAL "${${BUILD_FLAGS_VAR}}")
set(CXXFLAGS "${BUILD_FLAGS_VAL} ${CMAKE_CXX_FLAGS}")
set(HAVE_CMAKE 1 CACHE STRING "Whether built with CMake")

######################################################################
#
# Always use rpath to greatest extent.
# See: http://www.itk.org/Wiki/CMake_RPATH_handling
# Add -DCMAKE_SKIP_RPATH:BOOL=TRUE to prevent any rpath handling
#
######################################################################

# SET(CMAKE_SKIP_BUILD_RPATH FALSE)
# Use, i.e. don't skip the full RPATH for the build tree
# Not needed since we build static libs

# SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# When building, don't link with the install RPATH at build,
# but add later on when installing.
# Not needed since we build static libs

# SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# The RPATH to be used when installing
# Not needed since we build static libs

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Add the automatically determined parts of the RPATH that
# point to directories outside the build tree to the install RPATH

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

# z-lib only on windows
if (NOT WIN32)
  find_package(TxZ REQUIRED)
endif ()

# find HDF5
find_package(TxHdf5 REQUIRED)
if (HDF5_FOUND)
   include_directories(${Hdf5_INCLUDE_DIRS})
endif (HDF5_FOUND)

# find GSL
find_package(TxGsl REQUIRED)
if (GSL_FOUND)
   include_directories(${Gsl_INCLUDE_DIRS})
endif (GSL_FOUND)

# find FFTW3
find_package(TxFftw3)

# find Boost
find_package(TxBoost REQUIRED)
if (BOOST_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif (BOOST_FOUND)

# Math packages
find_package(TxPetsc REQUIRED)

IF(APPLE)
   # add the Accelerate framework on Mac OS X
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   MARK_AS_ADVANCED(ACCELERATE_LIB)
   SET(EXTRA_LIBS ${ACCELERATE_LIB})
else(APPLE)
  # CLapack is used on Unix and Windows
  # WARNING: On Unix should use system Lapack
  find_package(TxCLapackCMake REQUIRED)
  include_directories(${CLapackCMake_INCLUDE_DIRS})
endif(APPLE)

if (ENABLE_PARALLEL AND MPI_INCLUDE_PATH AND NOT TX_HAVE_MPI_COMPILER_WRAPPER)
  include_directories("${MPI_INCLUDE_PATH}")
endif ()

# Create configure file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in
   ${CMAKE_CURRENT_BINARY_DIR}/config.h)

######################################################################
#
# Set permissions before adding subdirectories
#
######################################################################

set(TX_GROUP_WRITE GROUP_WRITE)
set(TX_WORLD_FILE_PERMS WORLD_WRITE)
set(TX_WORLD_PROGRAM_PERMS WORLD_WRITE WORLD_EXECUTE)

install(FILES svninfo.txt DESTINATION include)

# add subdirectories
add_subdirectory(etc/loki/src)
add_subdirectory(etc/lua/src)
add_subdirectory(lib)
add_subdirectory(grids)
add_subdirectory(datastructs)
add_subdirectory(slvrs)
add_subdirectory(proto)
add_subdirectory(rte)
add_subdirectory(gkeyll)
add_subdirectory(unit)

# Add individual tests
ENABLE_TESTING()

# add tests
add_test(lcarray unit/lcarray)
add_test(lccolmajorindexer unit/lccolmajorindexer)
add_test(lcrowmajorindexer unit/lcrowmajorindexer)
add_test(lcfixedvector unit/lcfixedvector)
add_test(lcvector unit/lcvector)
add_test(lcmatrix unit/lcmatrix)
add_test(lccmdlineargs unit/lccmdlineargs)
add_test(lcregion unit/lcregion)
add_test(lcfield unit/lcfield)
add_test(lcrowmajorsequencer unit/lcrowmajorsequencer)
add_test(lccolmajorsequencer unit/lccolmajorsequencer)
add_test(lcdirsequencer unit/lcdirsequencer)
add_test(lcluatable unit/lcluatable)
add_test(lcluamodule unit/lcluamodule)
add_test(lclogger unit/lclogger)
add_test(lchdf5io unit/lchdf5io)
add_test(lcvec3 unit/lcvec3)
add_test(lcrectcartgrid unit/lcrectcartgrid)
add_test(lclincombiner unit/lclincombiner)
add_test(lcfieldptr unit/lcfieldptr)
add_test(lcrectcoordsys unit/lcrectcoordsys)
add_test(lcgeomprim unit/lcgeomprim)
add_test(lcstructgridfield unit/lcstructgridfield)
add_test(lcdecompregion unit/lcdecompregion)
add_test(lcmultiregion unit/lcmultiregion)
add_test(lcquadraturerule unit/lcquadraturerule)
add_test(lcdynvector unit/lcdynvector)

# special tests run via executables
if (ENABLE_PARALLEL)
  add_test(lcstructgridfield-p4 ${CMAKE_SOURCE_DIR}/unit/parunitexec.sh 4 ./unit/lcstructgridfield)
  add_test(lcstructgridfield-p2 ${CMAKE_SOURCE_DIR}/unit/parunitexec.sh 2 ./unit/lcstructgridfield)
endif (ENABLE_PARALLEL)

######################################################################
#
# Finalize
#
######################################################################

set(COMMON_INSTALL TRUE)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "LUCEE")
set(CONFIG_FILES)  # Clear out
include(${CMAKE_SOURCE_DIR}/CMake/TxFinalize.cmake)
