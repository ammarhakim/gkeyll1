# - FindSciBoost: Module to find include directories and libraries for
#   Boost. This module was originally developed to set the
#   SUPRA_SEARCH_PATH, so that the system path was not checked before
#   the user specified path, and included the stock FindBoost. This
#   changed after quite a few modifications, as it would still look at
#   system path if the libraries weren't found in the user specified
#   path.
#
#   Should be modified to include the stock FindBoost?
#
# This module can be included in CMake builds using find_package:
#   find_package(SciBoost REQUIRED signals filesystem system ...)
#
# The components list needs to contain actual names of boost libraries
# only: signals for libboost_signals, system for libboost_system, etc.
#
# This module will define the following variables:
#   BOOST_FOUND, HAVE_BOOST = True if Boost is found: the include directory was
#                             found and all the libraries specified were found.
#   Boost_INCLUDE_DIRS      = Location of Boost includes
#   Boost_LIBRARY_DIRS      = Location of Boost libraries
#   Boost_LIBRARIES         = Required libraries

######################################################################
#
# FindSciBoost: find includes and libraries for boost
#
# $Id: FindSciBoost.cmake 626 2014-09-10 16:04:55Z jrobcary $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

# Default: libraries have boost_ prepended.
set(BOOST_LIB_PREFIX boost_)
if (USE_CC4PY_LIBS)
  set(instdirs boost-cc4py boost-sersh)
# Shared windows boost has libboost_ prepended to the name
  set(BOOST_LIB_PREFIX boost_)
elseif (USE_SHARED_LIBS OR BUILD_SHARED_LIBS OR ENABLE_SHARED)
  set(instdirs boost-sersh)
  set(BOOST_LIB_PREFIX boost_)
elseif (BUILD_WITH_SHARED_RUNTIME)
  set(instdirs boost-sermd boost-sersh)
  if (WIN32)
    set(BOOST_LIB_PREFIX libboost_)
  endif ()
else ()
  message(STATUS "Setting boost up for static linking")
  set(Boost_USE_STATIC_RUNTIME OFF)
  set(Boost_USE_STATIC_LIBS ON)
  set(Boost_USE_MULTITHREADED ON)
  set(instdirs boost)
# Static cases Windows has libboost_ prepended to the name
  if (WIN32)
    set(BOOST_LIB_PREFIX libboost_)
  endif ()
endif ()
message(STATUS "instdirs = ${instdirs}.")

include(FindPackageHandleStandardArgs)
# Set names and dirs for finding boost
if (DEBUG_CMAKE)
  message(STATUS "SciBoost_FIND_COMPONENTS = ${SciBoost_FIND_COMPONENTS}.")
  message(STATUS "SciBoost_FIND_REQUIRED = ${SciBoost_FIND_REQUIRED}.")
endif ()
set(SciBoost_LIBRARY_LIST "")
foreach (COMPONENT ${SciBoost_FIND_COMPONENTS})
  set(SciBoost_LIBRARY_LIST ${SciBoost_LIBRARY_LIST} ${BOOST_LIB_PREFIX}${COMPONENT})
endforeach ()
if (DEBUG_CMAKE)
  message(STATUS "SciBoost_LIBRARY_LIST = ${SciBoost_LIBRARY_LIST}.")
endif ()

SciFindPackage(PACKAGE "Boost"
  INSTALL_DIRS ${instdirs}
  HEADERS "boost/thread.hpp"
  LIBRARIES "${SciBoost_LIBRARY_LIST}"
)
unset(SciBoost_LIBRARY_LIST CACHE)
if (DEBUG_CMAKE)
  message(STATUS "Boost_DLLS = ${Boost_DLLS}.")
endif ()
if (Boost_DLLS)
  message(STATUS "Correcting Boost DLL definition")
  set(Boost_DEFINITIONS -DBOOST_ALL_DYN_LINK)
endif ()

set(Boost_DEFINITIONS ${Boost_DEFINITIONS} -DBOOST_ALL_NO_LIB)  # Disable auto-linking
message(STATUS "Boost_DEFINITIONS = ${Boost_DEFINITIONS}.")

if (BOOST_FOUND AND NOT Boost_INCLUDE_DIRS)
  set(BOOST_FOUND FALSE)
  message(STATUS "Reversing Boost found as Boost_INCLUDE_DIRS is empty.")
endif ()

if (BOOST_FOUND)
  # message(STATUS "Found Boost")
  set(HAVE_BOOST 1 CACHE BOOL "Whether have the Boost library")
else ()
  message(STATUS "Did not find Boost.  Use -DBoost_ROOT_DIR to specify installation directory.")
  if (SciBoost_FIND_REQUIRED)
    message(FATAL_ERROR "Failed finding boost.")
  endif ()
endif ()

if ("${BOOST_LIB_PREFIX}" STREQUAL "libboost_")
  foreach (COMPONENT ${SciBoost_FIND_COMPONENTS})
    set(Boost_boost_${COMPONENT}_LIBRARY ${Boost_${BOOST_LIB_PREFIX}${COMPONENT}_LIBRARY})
  endforeach ()
endif ()



