# - FindSciSquish: Module to find executables and libraries
#   for Squish. A stock CMake FindSquish.cmake does exist but
#   is supported only for Squish version 3.
#
# This module can be included in CMake builds in find_package:
#   find_package(SciSquish REQUIRED)
#
# This module will define the following variables:
#  HAVE_SQUISH         = Whether have the Squish library
#  Squish_PROGRAMS = Location of Squish executables
#  Squish_LIBRARY_DIRS = Location of Squish libraries
#  Squish_LIBRARIES    = Required libraries

######################################################################
#
# SciFindSquish: find includes and libraries for Squish.
#
# $Id: FindSciSquish.cmake 414 2013-12-17 18:47:54Z techxdave $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################
set(SUPRA_SEARCH_PATH ${SUPRA_SEARCH_PATH})

SciFindPackage(PACKAGE "Squish"
              INSTALL_DIR "squish"
              FILES ".squish-3-license"
              PROGRAMS "squishrunner;squishserver;nchcompare"
              LIBRARIES "squishqtpre"
)

if (SQUISH_FOUND)
  message(STATUS "Found Squish")
  set(HAVE_SQUISH 1 CACHE BOOL "Whether have the SQUISH package")
else ()
  message(STATUS "Did not find Squish.  Use -DSQUISH_DIR to specify the installation directory.")
endif ()

