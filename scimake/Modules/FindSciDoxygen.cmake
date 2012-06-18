# - FindSciDoxygen: Module to find doxygen and setup apidocs target for
#   Doxygen.
#
# Module usage:
#   find_package(SciDoxygen ...)
#
# This module will define the following variables:
#  DOXYGEN_FOUND         = Whether Doxygen was found
#  DOXYGEN_EXECUTABLE    = Path to doxygen executable

######################################################################
#
# SciDoxygen: Find Doxygen and set up apidocs target
#
# $Id: FindSciDoxygen.cmake 1351 2012-06-01 11:04:59Z mdurant $
#
# Copyright 2011 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

message("")
message("--------- FindSciDoxygen looking for doxygen ---------")
find_package(Doxygen)

if (NOT DOXYGEN_FOUND)
  message(STATUS "Doxygen not found by default CMake module, falling back to SciFindPackage module.")
  SciFindPackage(PACKAGE Doxygen
    EXECUTABLES "doxygen"
  )
  if (DOXYGEN_FOUND)
    set(DOXYGEN_EXECUTABLE "${Doxygen_doxygen}")
  endif ()
endif ()

if (DOXYGEN_FOUND)
  message(STATUS "DOXYGEN_EXECUTABLE found.")
  message(STATUS "DOXYGEN_EXECUTABLE = ${DOXYGEN_EXECUTABLE}")
else ()
  message(STATUS "DOXYGEN_EXECUTABLE not found. API documentation cannot be built.")
  set(ENABLE_DEVELDOCS FALSE)
endif ()

