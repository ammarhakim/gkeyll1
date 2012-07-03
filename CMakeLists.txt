# Project information
project(gkeyll)

set(VERSION_MAJOR "0")
set(VERSION_MINOR "5")
set(VERSION_PATCH "0")
set(PROJECT_VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})
set(GKEYLL_VERSION "${PROJECT_VERSION}")

# Required version
cmake_minimum_required(VERSION 2.8.6)

######################################################################
#
# Usual initialization stuff
#
######################################################################

include(${CMAKE_SOURCE_DIR}/scimake/SciInit.cmake)
set(GKEYLL_REV "${PROJECT_REV}")
set(GKEYLL_URL ${PROJECT_URL})
set(GKEYLL_INSTALLDIR ${GKEYLL_INSTALLDIR})

######################################################################
#
# Set for inclusion in config.h
#
######################################################################

set(CXX ${CMAKE_CXX_COMPILER})
set(BUILD_FLAGS_VAR  CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
set(BUILD_FLAGS_VAL "${${BUILD_FLAGS_VAR}}")
set(CXXFLAGS "${BUILD_FLAGS_VAL} ${CMAKE_CXX_FLAGS}")
set(HAVE_CMAKE 1 CACHE STRING "Whether built with scimake")

######################################################################
#
# Always use rpath to greatest extent.
# See: http://www.itk.org/Wiki/CMake_RPATH_handling
# Add -DCMAKE_SKIP_RPATH:BOOL=TRUE to prevent any rpath handling
#
######################################################################

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# Add the automatically determined parts of the RPATH that
# point to directories outside the build tree to the install RPATH

######################################################################
#
# Find packages
#
######################################################################

# Prerequisites for finding packages
include(${CMAKE_SOURCE_DIR}/scimake/Modules/SciFindPackage.cmake)

# find txbase
find_package(SciTxBase REQUIRED)
if (TXBASE_FOUND)
   include_directories(${TxBase_INCLUDE_DIRS})
endif (TXBASE_FOUND)

# z-lib only on windows
if (NOT WIN32)
  find_package(SciZ REQUIRED)
endif ()

# find HDF5
find_package(SciHdf5 REQUIRED)
if (HDF5_FOUND)
   include_directories(${Hdf5_INCLUDE_DIRS})
endif (HDF5_FOUND)

# find GSL
find_package(SciGsl REQUIRED)
if (GSL_FOUND)
   include_directories(${Gsl_INCLUDE_DIRS})
endif (GSL_FOUND)

# find FFTW3
find_package(SciFftw3)

# find Boost
find_package(SciBoost REQUIRED)
if (BOOST_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif (BOOST_FOUND)

# find LUA
find_package(SciLua REQUIRED)
if (LUA_FOUND)
  include_directories(${Lua_INCLUDE_DIRS})
endif (LUA_FOUND)

# find LUABIND
find_package(SciLuabind REQUIRED)
if (LUABIND_FOUND)
  include_directories(${Luabind_INCLUDE_DIRS})
endif (LUABIND_FOUND)

# find Blitz
find_package(SciBlitz REQUIRED)
if (BLITZ_FOUND)
  include_directories(${Blitz_INCLUDE_DIRS})
endif (BLITZ_FOUND)

# find Eigen
find_package(SciEigen REQUIRED)
if (EIGEN_FOUND)
  include_directories(${Eigen_INCLUDE_DIRS})
endif (EIGEN_FOUND)

# Math packages
find_package(SciPetsc REQUIRED)

IF(APPLE)
   # add the Accelerate framework on Mac OS X
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   MARK_AS_ADVANCED(ACCELERATE_LIB)
   SET(EXTRA_LIBS ${ACCELERATE_LIB})
else(APPLE)
  # Lapack is used on Unix and Windows
  find_package(SciLapack REQUIRED)
  include_directories(${Lapack_INCLUDE_DIRS})
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

set(SCI_GROUP_WRITE GROUP_WRITE)
set(SCI_WORLD_FILE_PERMS WORLD_READ)
set(SCI_WORLD_PROGRAM_PERMS WORLD_READ WORLD_EXECUTE)

# add subdirectories
add_subdirectory(etc/loki/src)
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
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "GKEYLL Solver System")
set(CONFIG_FILES)  # Clear out
include(${CMAKE_SOURCE_DIR}/scimake/SciFinalize.cmake)
