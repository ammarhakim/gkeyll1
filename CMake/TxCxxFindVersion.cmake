######################################################################
#
# TxCxxFindVersion: Determine compiler version for any compiler
#
# $Id: TxCxxFindVersion.cmake 158 2010-11-16 18:21:53Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

IF(CMAKE_CXX_COMPILER MATCHES "icl")
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    OUTPUT_VARIABLE _cxx_VERSION
  )
  STRING(REGEX MATCH
    "w_cproc_p_[0-9][0-9]\\.[0-9]\\.[0-9][0-9][0-9]"
    _cxx_VERSION
    ${_cxx_VERSION}
  )
  STRING(REPLACE "w_cproc_p_" "" _cxx_VERSION ${_cxx_VERSION})
ELSEIF(CMAKE_CXX_COMPILER MATCHES "cl")
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    OUTPUT_VARIABLE _cxx_VERSION
  )
  STRING(REGEX MATCH
    "Version [0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+ for"
    _cxx_VERSION
    ${_cxx_VERSION}
  )
  STRING(REPLACE "Version " "" _cxx_VERSION ${_cxx_VERSION})
  STRING(REPLACE " for" "" _cxx_VERSION ${_cxx_VERSION})
ELSEIF(${CMAKE_CXX_COMPILER_ID} STREQUAL PathScale)
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    ARGS --version
    OUTPUT_VARIABLE _cxx_VERSION
  )
  # ARGS -v # This used to work above?
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
  STRING(REGEX MATCH
    "Version [0-9]+\\.[0-9]"
    _cxx_VERSION
    ${_cxx_VERSION}
  )
  STRING(REPLACE "Version " "" _cxx_VERSION ${_cxx_VERSION})
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
ELSEIF(${CMAKE_CXX_COMPILER_ID} STREQUAL PGI)
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    ARGS -V
    OUTPUT_VARIABLE _cxx_VERSION
  )
  STRING(REGEX MATCH
    "pgCC [0-9]+\\.[0-9]+-[0-9]+"
    _cxx_VERSION
    ${_cxx_VERSION}
  )
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
  STRING(REPLACE "pgCC " "" _cxx_VERSION ${_cxx_VERSION})
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
ELSEIF(${CMAKE_CXX_COMPILER_ID} STREQUAL XL)
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    ARGS -qversion
    OUTPUT_VARIABLE _cxx_VERSION
  )
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
  STRING(REGEX MATCH
    "Version: .*"
    _cxx_VERSION
    ${_cxx_VERSION}
  )
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
  STRING(REPLACE "Version: " "" _cxx_VERSION ${_cxx_VERSION})
  # MESSAGE("_cxx_VERSION = ${_cxx_VERSION}.")
ELSE()
  INCLUDE(${CMAKE_ROOT}/Modules/FindBoost.cmake)
  _Boost_COMPILER_DUMPVERSION(_cxx_VERSION)
ENDIF(CMAKE_CXX_COMPILER MATCHES "icl")

set(CXX_VERSION ${_cxx_VERSION})

