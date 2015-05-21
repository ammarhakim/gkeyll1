# - FindSciGsl: Module to find include directories and
#   libraries for Gsl.
#
# Module usage:
#   find_package(SciCgns ...)
#
# This module will define the following variables:
#  HAVE_GSL, GSL_FOUND = Whether libraries and includes are found
#  Cgns_INCLUDE_DIRS       = Location of Gsl includes
#  Cgns_LIBRARY_DIRS       = Location of Gsl libraries
#  Cgns_LIBRARIES          = Required libraries

######################################################################
#
# FindCgns: find includes and libraries for GSL
#
# $Id: FindSciCgns.cmake 792 2015-04-17 14:07:44Z jrobcary $
#
# Copyright 2010-2015, Tech-X Corporation, Boulder, CO.
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#
######################################################################

SciFindPackage(PACKAGE "Cgns"
              INSTALL_DIR "cgns"
              HEADERS "cgnslib.h"
              LIBRARIES "cgns"
              INCLUDE_SUBDIRS "include"
              LIBRARY_SUBDIRS "lib"
             )

if (GSL_FOUND)
  message(STATUS "Found Cgns")
  set(HAVE_CGNS 1 CACHE BOOL "Whether have Cgns")
else ()
  message(STATUS "Did not find Cgns.  Use -DCGNS_DIR to specify the installation directory.")
  if (SciCgns_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
