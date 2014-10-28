# - FindSciTranspbase: Module to find include directories and libraries
#   for Transpbase. This module was implemented as there is no stock
#   CMake module for Transpbase.
#
# This module can be included in CMake builds in find_package:
#   find_package(SciTranspbase REQUIRED)
#
# This module will define the following variables:
#  HAVE_TRANSPBASE         = Whether have the Transpbase library
#  Transpbase_INCLUDE_DIRS = Location of Transpbase includes
#  Transpbase_LIBRARY_DIRS = Location of Transpbase libraries
#  Transpbase_LIBRARIES    = Required libraries
#  Transpbase_STLIBS       = Location of Transpbase static library

######################################################################
#
# FindSciTranspbase: find includes and libraries for Transpbase
#
# $Id: FindSciTranspbase.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

if (ENABLE_PARALLEL)
  set(instdir transpbase-par)
else ()
  set(instdir transpbase)
endif ()
SciFindPackage(PACKAGE "Transpbase"
  INSTALL_DIR "${instdir}"
  HEADERS "transperror.h"
  LIBRARIES "transpbase"
  LIBRARY_SUBDIRS "lib/${Fortran_COMP_LIB_SUBDIR};lib"
)

if (TRANSPBASE_FOUND)
  # message(STATUS "Found TRANSPBASE")
  set(HAVE_TRANSPBASE 1 CACHE BOOL "Whether have the TRANSPBASE libraries")
endif ()

