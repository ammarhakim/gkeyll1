######################################################################
#
# TxFindPackage: find includes and libraries of a package
#
# $Id: TxFindPackage.cmake 154 2010-10-28 14:28:45Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

# TxGetDir
#
# Args:
# canddirvar: the variable holding the directory candidate
# realdirvar: the variable holding the directory after following shortcuts, if needed
FUNCTION(TxGetDir canddirvar realdirvar)
  # MESSAGE("realdirvar = ${realdirvar}.")
  SET(canddir ${${canddirvar}})
  SET(${realdirvar})
  IF(canddir)
    IF(DEBUG_CMAKE)
      MESSAGE("Finding the directory for ${canddir}.")
    ENDIF()
    IF(WIN32)
      # MESSAGE("Looking for ${canddir}.lnk")
      IF(EXISTS ${canddir}.lnk)
        # MESSAGE("${canddir}.lnk exists.")
        EXEC_PROGRAM(readshortcut ARGS ${canddir}.lnk
          OUTPUT_VARIABLE idir
          RETURN_VALUE ires)
        IF(ires)
          MESSAGE(FATAL_ERROR "readshortcut error on ${canddir}.lnk.")
        ENDIF(ires)
        EXEC_PROGRAM(cygpath ARGS -am ${idir}
          OUTPUT_VARIABLE rd)
      ELSE(EXISTS ${canddir}.lnk)
        SET(rd ${canddir})
      ENDIF(EXISTS ${canddir}.lnk)
    ELSE(WIN32)
      SET(rd ${canddir})
    ENDIF(WIN32)
    # MESSAGE("Found ${rd}.")
  ENDIF(canddir)
  SET(${realdirvar} ${rd} PARENT_SCOPE)
ENDFUNCTION(TxGetDir canddir realdir)

# TxFindPackage
#
# Args:
# txpkg: the prefix for the defined variables
# txinst: the name for the installation directory.  Defaults
#   to lower cased txpkgname
# txhdrs: the header files to look for.
# txlibs: the libraries
# txmods: Fortran module files
# txincsd: include subdirectories
# txlibsd: library subdirectories
#
# NOTE: One cannot reset calling variables
MACRO(TxFindPackage txpkg txinst
  txhdrs txlibs txmods
  txincsd txlibsd)

  IF(DEBUG_CMAKE)
    MESSAGE(STATUS "${TxFindPackage}
    ${txpkg}
    ${txinst}
    ${txhdrs}
    ${txlibs}
    ${txmods}
    ${txincsd}
    ${txlibsd}")
  ENDIF()
# Construct variable name
  STRING(REGEX REPLACE "[./-]" "_" txpkgreg ${txpkg})
# txpkgreg is the regularized package name
  STRING(TOUPPER ${txpkgreg} txpkguc)
  IF(DEBUG_CMAKE)
    MESSAGE(STATUS "ENABLE_${txpkguc} = ${ENABLE_${txpkguc}}")
  ENDIF()
  IF(NOT DEFINED ENABLE_${txpkguc})
    SET(ENABLE_${txpkguc} TRUE)
  ENDIF()
  IF(NOT ENABLE_${txpkguc})
    MESSAGE(STATUS "Disabling ${txpkgreg}.")
    SET(${txpkguc}_FOUND FALSE)
    IF(DEBUG_CMAKE)
      MESSAGE(STATUS "${txpkguc}_FOUND set to FALSE.")
    ENDIF()
    RETURN()
  ENDIF()
  STRING(TOLOWER ${txpkgreg} txpkglc)
  # MESSAGE("txpkgreg = ${txpkgreg}")
  # MESSAGE("txpkguc = ${txpkguc}")
  # MESSAGE("txpkglc = ${txpkglc}")
  # MESSAGE("txinst = ${txinst}")
  STRING(LENGTH "${txinst}" txinstlen)
  #MESSAGE("txinstlen = ${txinstlen}")
  # IF("${txinst}")  # This does not work
  # IF(txinst)       # This does not work
  IF(${txinstlen})   # This works.  Why so difficult?
    # MESSAGE("txinst is defined.")
    SET(txpkginst "${txinst}")
  ELSE()
    # MESSAGE("txinst is not defined.")
    SET(txpkginst ${txpkglc})
  ENDIF()

# Find the possible directories
# Find the translated supra search path
  SET(txpath)
  # SET(txdir ${${txpkgreg}_DIR})
  SET(txdir ${${txpkguc}_DIR})
  # MESSAGE("Calling TxGetDir on ${txdir}.")
  TxGetDir(txdir txpath)
  FOREACH(instdir ${txpkginst})
    FOREACH(spdir ${SUPRA_SEARCH_PATH})
      SET(idir ${spdir}/${instdir})
      TxGetDir(idir txdir)
      SET(txpath ${txpath} ${txdir})
    ENDFOREACH(spdir ${SUPRA_SEARCH_PATH})
  ENDFOREACH()
  LIST(LENGTH txpath txpathlen)
  IF(DEBUG_CMAKE)
    IF (NOT txpathlen)
      MESSAGE(FATAL_ERROR "txpath is empty.")
    ELSE()
      MESSAGE(STATUS "txpath = ${txpath}")
    ENDIF()
  ENDIF()

# Finding none is okay as in finding zlib.
  STRING(LENGTH "${txhdrs}" txhdrslen)
  IF(${txhdrslen})
# Find include subdirs
    STRING(LENGTH "${txincsd}" txlen)
    IF(${txlen})
      SET(txincsubdirs ${txincsd})
    ELSE()
      SET(txincsubdirs include)
    ENDIF()
    IF(DEBUG_CMAKE)
      MESSAGE("txincsubdirs = ${txincsubdirs}.")
    ENDIF()
    IF(DEBUG_CMAKE)
      MESSAGE("Looking for headers and libs under ${txpath} with subdirs, ${txincsubdirs},")
    ENDIF()
    SET(${txpkgreg}_INCLUDE_DIR)
    FOREACH(txhdr ${txhdrs})
# First look in specified path
      UNSET(${txpkgreg}_INCLUDE_DIR)
      FIND_PATH(${txpkgreg}_INCLUDE_DIR
        ${txhdr}
        PATHS ${txpath}
        PATH_SUFFIXES "${txincsubdirs}"
        NO_DEFAULT_PATH)
      IF(NOT ${txpkgreg}_INCLUDE_DIR)
        FIND_PATH(${txpkgreg}_INCLUDE_DIR
          ${txhdr}
          PATHS ${txpath}
          PATH_SUFFIXES "${txincsubdirs}")
      ENDIF()
      STRING(REGEX REPLACE "[./-]" "_" txhdrvar ${txhdr})
      # STRING(TOUPPER ${txhdrvar} txhdrvar)
      SET(txhdrvar ${txpkgreg}_${txhdrvar})
      SET(${txhdrvar} "${${txpkgreg}_INCLUDE_DIR}/${txhdr}")
      # MESSAGE("${txhdrvar} = ${${txhdrvar}}")
      SET(${txpkgreg}_INCLUDES ${${txpkgreg}_INCLUDES}
	${${txpkgreg}_INCLUDE_DIR})
    ENDFOREACH(txhdr ${txhdrs})
    LIST(REMOVE_DUPLICATES ${txpkgreg}_INCLUDES)
    # MESSAGE(STATUS "${txpkgreg}_INCLUDES = ${${txpkgreg}_INCLUDES}")
    LIST(LENGTH ${txpkgreg}_INCLUDES txinclistlen)
    IF(NOT ${txinclistlen})
      MESSAGE(WARNING "None of ${txhdrs} found.  Define ${txpkguc}_DIR to find them.")
    ENDIF(NOT ${txinclistlen})
  ENDIF(${txhdrslen})

# Find libraries.
# Finding none is fatal, but finding a subset is okay as allows us
# to look for a maximum set, as needed for Trilinos.
  # MESSAGE("txlibs = ${txlibs}")
  STRING(LENGTH "${txlibs}" txlibslen)
  # MESSAGE("txlibsd = ${txlibsd}")
  IF(${txlibslen})
# Add in extrax subdirs
    STRING(LENGTH "${txlibsd}" txlen)
    IF(${txlen})
      SET(txlibsubdirs ${txlibsd})
    ELSE(${txlen})
      SET(txlibsubdirs lib)
    ENDIF(${txlen})
# If library_dirs specified, add that to the front of the path
    IF(${txpkgreg}_LIBRARY_DIRS)
      SET(txlibsubdirs . ${txlibsubdirs})
      SET(txpath ${${txpkgreg}_LIBRARY_DIRS} ${txpath})
    ENDIF(${txpkgreg}_LIBRARY_DIRS)
# Look for this library
    SET(${txpkgreg}_LIBRARY_DIR)
    FOREACH(txlib ${txlibs})
      STRING(REGEX REPLACE "[./-]" "_" txlibvar ${txlib})
      SET(txlibvar ${txpkgreg}_${txlibvar}_LIBRARY)
# First look in defined path
      UNSET(${txlibvar})
      IF(DEBUG_CMAKE)
        MESSAGE("Looking for ${txlib} in ${txpath} with subdirs, ${txlibsubdirs}.")
      ENDIF()
      FIND_LIBRARY(${txlibvar}
        ${txlib}
        PATHS ${txpath}
        PATH_SUFFIXES ${txlibsubdirs}
        NO_DEFAULT_PATH)
      IF(DEBUG_CMAKE)
        MESSAGE("After initial search, ${txlibvar} = ${${txlibvar}}.")
      ENDIF()
      IF(NOT ${txlibvar})
        FIND_LIBRARY(${txlibvar}
          ${txlib}
          PATHS ${txpath}
          PATH_SUFFIXES ${txlibsubdirs})
      ENDIF()
      IF(DEBUG_CMAKE)
        MESSAGE("After second search, ${txlibvar} = ${${txlibvar}}.")
      ENDIF()
      IF(${txlibvar})
        SET(${txpkgreg}_LIBRARIES ${${txpkgreg}_LIBRARIES} ${${txlibvar}})
      ENDIF()
    ENDFOREACH(txlib ${txlibs})
    LIST(LENGTH ${txpkgreg}_LIBRARIES txliblistlen)
    IF(NOT ${txliblistlen})
      MESSAGE(STATUS "None of the libraries, ${txlibs}, found.  Define ${txpkguc}_DIR to find them.")
    ELSE(NOT ${txliblistlen})
      LIST(REMOVE_DUPLICATES ${txpkgreg}_LIBRARIES)
      # MESSAGE(STATUS "${txpkgreg}_LIBRARIES = ${${txpkgreg}_LIBRARIES}")
    ENDIF(NOT ${txliblistlen})
  ELSE(${txlibslen})
    SET(txliblistlen 0)
  ENDIF(${txlibslen})

# If libraries were requested, and any were found, or
# if any headers were found, then we have this.
  # IF(${txlibslen} AND ${txliblistlen} OR ${txhdrslen})
  IF(${${txpkgreg}_INCLUDES} MATCHES ${txpkgreg}_INCLUDE_DIR-NOTFOUND)
    MESSAGE("${txpkguc} not found.")
  ELSE()
    SET(${txpkguc}_FOUND TRUE)
    IF(DEBUG_CMAKE)
      MESSAGE("${txpkguc}_FOUND set to TRUE")
    ENDIF()
  ENDIF()

# Signal the results
# Per http://www.cmake.org/Wiki/CMake:How_To_Find_Installed_Software,
# The convention is to capitalize the _FOUND variable.
  IF(${txpkguc}_FOUND)
    IF(NOT ${txpkgreg}_FIND_QUIETLY)
      MESSAGE(STATUS "Found ${txpkg}")
      MESSAGE(STATUS "${txpkgreg}_INCLUDES = ${${txpkgreg}_INCLUDES}")
      MESSAGE(STATUS "${txpkgreg}_LIBRARIES = ${${txpkgreg}_LIBRARIES}")
    ENDIF(NOT ${txpkgreg}_FIND_QUIETLY)
  ELSE(${txpkguc}_FOUND)
    IF(${txpkgreg}_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find ${txpkg}")
    ENDIF(${txpkgreg}_FIND_REQUIRED)
  ENDIF(${txpkguc}_FOUND)

ENDMACRO(TxFindPackage)

