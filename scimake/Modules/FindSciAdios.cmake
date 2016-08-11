# - FindSciAdios: Module to find include directories and libraries for
#   Adios. This module was implemented as there is no stock CMake
#   module for Adios.
#
# This module can be included in CMake builds in find_package:
#   find_package(SciAdios REQUIRED)
#
# This module will define the following variables:
#  HAVE_ADIOS         = Whether have the Adios library
#  Adios_INCLUDE_DIRS = Location of Adios includes
#  Adios_LIBRARY_DIRS = Location of Adios libraries
#  Adios_LIBRARIES    = Required libraries

######################################################################
#
# SciFindAdios: find includes and libraries for Adios
#
# $Id: FindSciAdios.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

SciFindPackage(PACKAGE "Adios"
              INSTALL_DIR "adios"
              HEADERS "adios.h"
              LIBRARIES "adios"
              MODULES "adios"
              )

if (ADIOS_FOUND)
  message(STATUS "Found Adios")
  set(HAVE_ADIOS 1 CACHE BOOL "Whether have the ADIOS library")
else ()
  message(STATUS "Did not find Adios.  Use -DADIOS_DIR to specify the installation directory.")
  if (SciAdios_FIND_REQUIRED)
    message(FATAL_ERROR "Failed")
  endif ()
endif ()

