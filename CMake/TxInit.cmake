######################################################################
#
# TxInit: Do the startup stuff for any package
#
# $Id: TxInit.cmake 183 2010-12-29 20:12:40Z yjchoi $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

######################################################################
#
# Set up standard paths
#
######################################################################

# MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
# MESSAGE("CMAKE_ROOT = ${CMAKE_ROOT}")
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules ${CMAKE_ROOT}/Modules ${CMAKE_MODULE_PATH})
IF(DEBUG_CMAKE)
  MESSAGE("CMAKE_MODULE_PATH = ${CMAKE_MODULE_PATH}")
ENDIF()
IF(SUPRA_SEARCH_PATH)
  SET(SUPRA_SEARCH_PATH ${SUPRA_SEARCH_PATH})
  MESSAGE("SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")
ELSE(SUPRA_SEARCH_PATH)
  IF(WIN32)
    SET(SUPRA_SEARCH_PATH /cygwin $ENV{HOME}/software /cygwin/opt /cygwin/internal /cygwin/volatile /cygwin/contrib)
  ELSE(WIN32)
    SET(SUPRA_SEARCH_PATH $ENV{HOME}/software /opt /internal /volatile /contrib /usr/local /usr/lib /usr/bin)
    MESSAGE("SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")	
  ENDIF(WIN32)
ENDIF(SUPRA_SEARCH_PATH)
IF(DEBUG_CMAKE)
  MESSAGE("SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")
ENDIF()

######################################################################
#
# Get revisions
#
######################################################################

# MESSAGE("SVN_BIN = ${SVN_BIN}")
IF(SVN_BINDIR)
  SET(SVN_PATH ${SVN_BINDIR} ${PATH})
  MESSAGE("Looking for svn in ${SVN_PATH}")
  FIND_PROGRAM(SVN_BIN NAME svn PATHS ${SVN_PATH}
    DOC "subversion command line client" NO_DEFAULT_PATH)
  FIND_PROGRAM(SVNVERSION_BIN NAME svnversion PATHS ${SVN_PATH}
    DOC "subversion version command line client" NO_DEFAULT_PATH)
ENDIF()
# If not found, search in system paths
IF(NOT SVN_BIN)
  FIND_PROGRAM(SVN_BIN NAME svn PATHS ${SVN_PATH}
    DOC "subversion command line client")
  FIND_PROGRAM(SVNVERSION_BIN NAME svnversion PATHS ${SVN_PATH}
    DOC "subversion version command line client")
ENDIF(NOT SVN_BIN)
IF(DEBUG_CMAKE)
  MESSAGE(STATUS "SVN_BIN is ${SVN_BIN}")
  MESSAGE(STATUS "SVNVERSION_BIN is ${SVNVERSION_BIN}")
ENDIF()
IF(SVN_BIN AND SVNVERSION_BIN)
  INCLUDE(${PROJECT_SOURCE_DIR}/CMake/tx_svn_info.cmake)
  Subversion_GET_VERSION(${PROJECT_SOURCE_DIR} PROJECT_REV PROJECT_URL)
ENDIF(SVN_BIN AND SVNVERSION_BIN)
MESSAGE(STATUS "PROJECT_REV is ${PROJECT_REV}")
MESSAGE(STATUS "PROJECT_URL is ${PROJECT_URL}")
IF(SVN_BIN AND SVNVERSION_BIN)
Subversion_GET_VERSION(${PROJECT_SOURCE_DIR}/CMake CMAKEDIR_REV CMAKEDIR_URL)
ENDIF(SVN_BIN AND SVNVERSION_BIN)
MESSAGE(STATUS "CMAKEDIR_REV is ${CMAKEDIR_REV}")
MESSAGE(STATUS "CMAKEDIR_URL is ${CMAKEDIR_URL}")

######################################################################
#
# C, CXX Checks
#
######################################################################

IF(WIN32)
  ADD_DEFINITIONS(/DHAVE_CONFIG_H)
ELSE()
  ADD_DEFINITIONS(-DHAVE_CONFIG_H)
ENDIF()
INCLUDE(${PROJECT_SOURCE_DIR}/CMake/TxCChecks.cmake)
INCLUDE(${PROJECT_SOURCE_DIR}/CMake/TxCxxChecks.cmake)

######################################################################
#
# Fix shared flags on windows
#
######################################################################

INCLUDE(${PROJECT_SOURCE_DIR}/CMake/TxWinFlags.cmake)

######################################################################
#
# Look for MPI
#
######################################################################

INCLUDE(${PROJECT_SOURCE_DIR}/CMake/TxFindMpi.cmake)

######################################################################
#
# Load Find Package
#
######################################################################

INCLUDE(${PROJECT_SOURCE_DIR}/CMake/Modules/TxFindPackage.cmake)

######################################################################
#
# Host information
#
######################################################################

site_name(HOSTNAME)
MESSAGE(STATUS "hostname is ${HOSTNAME}")

