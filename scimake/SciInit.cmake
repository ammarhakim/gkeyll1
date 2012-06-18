######################################################################
#
# SciInit: Do the startup stuff for any package
#
# $Id: SciInit.cmake 1330 2012-04-23 18:36:31Z cary $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

#####################################################################
#
# Pull in useful macros
#
#####################################################################

if (NOT DEFINED SCICMAKE_DIR)
  set(SCICMAKE_DIR ${PROJECT_SOURCE_DIR}/scimake)
endif ()
include(${SCICMAKE_DIR}/SciFuncsMacros.cmake)
include(${SCICMAKE_DIR}/SciGetDepsFromInstall.cmake)

#####################################################################
#
# Clean out config.summary
#
#####################################################################

set(CONFIG_SUMMARY ${PROJECT_BINARY_DIR}/config.summary)
file(REMOVE ${CONFIG_SUMMARY})
SciPrintString("CONFIGURING ${CMAKE_PROJECT_NAME} with scimake in "
 "${PROJECT_BINARY_DIR}.")

#####################################################################
#
# Set some vars to upper case for case-free comparisons
#
#####################################################################

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UC)

#####################################################################
#
# Set group write installation perm
# May need to be different on different machines
#
#####################################################################

if (NOT DEFINED SCI_GROUP_WRITE)
  set(SCI_GROUP_WRITE ""
      CACHE STRING "Value of group write permissions")
endif ()
if (NOT DEFINED SCI_WORLD_PERMS)
  set(SCI_WORLD_PERMS ""
      CACHE STRING "Value of world permissions")
endif ()

#####################################################################
#
# Set OS-specific flags
#
#####################################################################

message(STATUS "[SciInit.cmake] CMAKE_SYSTEM_NAME is ${CMAKE_SYSTEM_NAME}")
if ("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  set(MACX TRUE CACHE BOOL "True if compiled on Mac OS X")
  message(STATUS "[SciInit.cmake] Compiling on MAC")
elseif ("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  set(LINUX TRUE CACHE BOOL "True if compiled on Linux")
  message(STATUS "[SciInit.cmake] Compiling on LINUX")
elseif (WIN32)
  set(WINDOWS TRUE CACHE BOOL "True if compiled on Windows")
  message(STATUS "[SciInit.cmake] Compiling on WINDOWS")
else ()
  message(FATAL_ERROR "[SciInit.cmake] Unrecognized OS!")
endif ()
SciPrintString("CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}.")

######################################################################
#
# Set up standard paths
#
######################################################################

# message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}.")
# message(STATUS "CMAKE_ROOT = ${CMAKE_ROOT}.")
# if (NOT DEFINED CMAKE_MODULE_PATH)
#  set(CMAKE_MODULE_PATH ${CMAKE_ROOT}/Modules)
# endif ()
# message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}.")
cmake_policy(SET CMP0017 OLD) # Use our modules over theirs
set(CMAKE_MODULE_PATH
  ${SCICMAKE_DIR}/Modules
)

set(TXCMAKE_DIR ${PROJECT_SOURCE_DIR}/txcmake)
if (EXISTS ${TXCMAKE_DIR})
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TXCMAKE_DIR}/Modules)
endif ()

if (DEBUG_CMAKE)
  message(STATUS "CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
endif ()

# message("SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")
if (SUPRA_SEARCH_PATH)
  set(SUPRA_SEARCH_PATH "${SUPRA_SEARCH_PATH}")
else ()
  if (WIN32)
# According to JRC the following must be turned off to compile under
# Windows(AP)
    # set(SUPRA_SEARCH_PATH $ENV{HOME}/software /winsame/internal /winsame/volatile /winsame/contrib /opt /usr/local)
  else ()
# JRC: SUPRA_SEARCH_PATH should include only top directory
# also, system paths should not be needed due to cmake's system search
    set(SUPRA_SEARCH_PATH $ENV{HOME}/software /internal /volatile /contrib /opt /usr/local)
  endif ()
endif ()
SciPrintString("SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")

######################################################################
#
# Get system description
#
######################################################################

find_program(HOSTNAME_CMD NAMES hostname)
exec_program(${HOSTNAME_CMD} ARGS OUTPUT_VARIABLE HOSTNAME)
SciPrintString("scimake running on ${HOSTNAME}")
string(REGEX REPLACE "\\..*$" "" UQHOSTNAME "${HOSTNAME}")
SciPrintString("UQHOSTNAME = ${UQHOSTNAME}")

find_program(UNAME NAMES uname)
macro(getuname name flag)
  exec_program("${UNAME}" ARGS "${flag}" OUTPUT_VARIABLE "${name}")
endmacro(getuname)

getuname(osname -s)
getuname(osrel  -r)
getuname(cpu    -m)
set(HOSTTYPE "${osname}-${cpu}")
SciPrintString("HOSTTYPE = ${HOSTTYPE}")
site_name(HOSTNAME)
SciPrintString("hostname is ${HOSTNAME}")

######################################################################
#
# Other useful provenance
#
######################################################################

# Calling it "ENV Var" + "NAME" to avoid conflicts elsewhere
set(HOMENAME $ENV{HOME})
set(USERNAME $ENV{USER})
set(SCRATCHNAME $ENV{SCRATCH})

######################################################################
#
# Get revisions
#
######################################################################

include(${SCICMAKE_DIR}/SciSvnInfo.cmake)

######################################################################
#
# config.h
#
######################################################################

if (NOT NO_CONFIG_H)
  if (WIN32)
    add_definitions(/DHAVE_CONFIG_H)
  else ()
    add_definitions(-DHAVE_CONFIG_H)
  endif ()
endif ()

######################################################################
#
# Fix shared flags on windows
#
######################################################################

include(${SCICMAKE_DIR}/SciWinFlags.cmake)

######################################################################
#
# C, CXX, Fortran Checks
#
######################################################################

include(${SCICMAKE_DIR}/SciCChecks.cmake)
if (NOT NOCXX)
  include(${SCICMAKE_DIR}/SciCxxChecks.cmake)
endif ()
if (NOT NOFORTRAN)
  message("")
  message("--------- SciFortranChecking ---------")
# Enable Fortran to all those variables
  enable_language(Fortran)
  include(${CMAKE_ROOT}/Modules/CMakeDetermineFortranCompiler.cmake)
  include(${SCICMAKE_DIR}/SciFortranChecks.cmake)
else ()
  message(STATUS "No Fortran, so no implicit fortran link libraries known.")
endif ()

######################################################################
#
# Need to add library paths from compiler for rpath
#
######################################################################


if ("${CMAKE_SYSTEM_NAME}" STREQUAL Linux)
  message("")
  message("--------- Checking for rpath ---------")
  # message("is Linux")
  if ("${CMAKE_C_COMPILER_ID}" STREQUAL GNU)
    # message("is GNU")
    execute_process(COMMAND ${CMAKE_C_COMPILER} -print-file-name=libstdc++.so
      OUTPUT_VARIABLE libcxx)
    message("libcxx is ${libcxx}")
    if (${libcxx} MATCHES "^/")
      get_filename_component(CXX_LIBDIR ${libcxx}/.. REALPATH)
      message(STATUS "libstdc++ is in ${CXX_LIBDIR}.")
# Add to build rpath
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,${CXX_LIBDIR}")
# For the installation rpath, the below works
      # set(CMAKE_INSTALL_RPATH ${CXX_LIBDIR})
# but it puts this in front on the installation rpath
# Since we want this after the origin stuff we instead add it to
# particular target properties as shown below in the particular dirs.
      # set_target_properties(mytarget
        # PROPERTIES INSTALL_RPATH "\$ORIGIN:\$ORIGIN/../lib:${CXX_LIBDIR}")
    endif ()
  endif ()
endif()

######################################################################
#
# For MinGW set libraries to Windows style
#
######################################################################

if (USING_MINGW)
  message("")
  message("--------- Setting MinGW library prefix and suffix to windows style  ---------")
  set(CMAKE_STATIC_LIBRARY_PREFIX "")
  set(CMAKE_STATIC_LIBRARY_SUFFIX .lib)
  message(STATUS "CMAKE_STATIC_LIBRARY_PREFIX = ${CMAKE_STATIC_LIBRARY_PREFIX}.")
  message(STATUS "CMAKE_STATIC_LIBRARY_SUFFIX = ${CMAKE_STATIC_LIBRARY_SUFFIX}.")
endif ()

######################################################################
#
# Load Find Package
#
######################################################################

include(${SCICMAKE_DIR}/Modules/SciFindPackage.cmake)

######################################################################
#
# Look for MPI
#
######################################################################

option(ENABLE_PARALLEL "Enable parallel build" OFF)
message("")
if (ENABLE_PARALLEL)
  message(STATUS "ENABLE_PARALLEL requested.  Will search for MPI.")
elseif (INSTALL_PARALLEL)
  message(STATUS "INSTALL_PARALLEL requested.  Will search for MPI.")
else ()
  message(STATUS "Not searching for MPI.")
endif ()
if (ENABLE_PARALLEL OR INSTALL_PARALLEL)
  find_package(SciMpi REQUIRED)
endif ()

######################################################################
#
# Host information
#
######################################################################


######################################################################
#
# Mac OS X make dylibs install with rpath
#
######################################################################

if (APPLE)
  set (CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX})
endif ()

