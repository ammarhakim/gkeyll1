# - FindSciMxml: Module to find include directories and
#   libraries for Mxml.
#
# Module usage:
#   find_package(SciMxml ...)
#
# This module will define the following variables:
#  HAVE_MXML, MXML_FOUND = Whether libraries and includes are found
#  Mxml_INCLUDE_DIRS       = Location of Mxml includes
#  Mxml_LIBRARY_DIRS       = Location of Mxml libraries
#  Mxml_LIBRARIES          = Required libraries

######################################################################
#
# SciFindMxml: find includes and libraries for txbase
#
# $Id: FindSciMxml.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

# mxml is only built in serial
SciFindPackage(PACKAGE "Mxml"
              INSTALL_DIR "mxml"
              HEADERS "mxml.h"
              INCLUDES "mxml"
              )

if (MXML_FOUND)
  message(STATUS "Found Mxml")
  set(HAVE_MXML 1 CACHE BOOL "Whether have the MXML library")
else ()
  message(STATUS "Did not find Mxml.  Use -DMXML_DIR to specify the installation directory.")
  if (SciMxml_FIND_REQUIRED)
    message(FATAL_ERROR "Failing")
  endif ()
endif ()


