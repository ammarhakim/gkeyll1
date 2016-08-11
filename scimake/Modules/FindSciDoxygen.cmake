# - FindSciDoxygen: Module to find doxygen and setup apidocs target for
#   Doxygen.
#
# Module usage:
#   find_package(SciDoxygen ...)
#
# This module will define the following variables:
#  DOXYGEN_FOUND         = Whether Doxygen was found
#  DOXYGEN_PROGRAM    = Path to doxygen executable

######################################################################
#
# SciDoxygen: Find Doxygen and set up apidocs target
#
# $Id: FindSciDoxygen.cmake 407 2013-12-14 19:04:56Z jrobcary $
#
# Copyright 2011 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

message("")
message("--------- FindSciDoxygen looking for doxygen ---------")
find_package(Doxygen)

if (NOT DOXYGEN_FOUND)
  message(STATUS "Doxygen not found by default CMake module, falling back to SciFindPackage module.")
  SciFindPackage(PACKAGE Doxygen
    PROGRAMS "doxygen"
  )
  if (DOXYGEN_FOUND)
    set(DOXYGEN_PROGRAM "${Doxygen_doxygen}")
  endif ()
endif ()

# Maintain backward compatibility
if (NOT DOXYGEN_PROGRAM)
  set(DOXYGEN_PROGRAM ${DOXYGEN_EXECUTABLE})
elseif (NOT DOXYGEN_EXECUTABLE)
  set(DOXYGEN_EXECUTABLE ${DOXYGEN_PROGRAM})
endif ()

if (DOXYGEN_FOUND)
  message(STATUS "DOXYGEN_PROGRAM found.")
  message(STATUS "DOXYGEN_PROGRAM = ${DOXYGEN_PROGRAM}")
  message(STATUS "DOXYGEN_EXECUTABLE = ${DOXYGEN_EXECUTABLE}")
else ()
  message(STATUS "DOXYGEN_PROGRAM not found. API documentation cannot be built.")
  set(ENABLE_DEVELDOCS FALSE)
endif ()

