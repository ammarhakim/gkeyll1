# - FindSciBlitz: Module to find include directories and
#   libraries for Blitz.
#
# Module usage:
#   find_package(SciBlitz ...)
#
# This module will define the following variables:
#  HAVE_BLITZ, BLITZ_FOUND = Whether libraries and includes are found
#  Blitz_INCLUDE_DIRS       = Location of Blitz includes
#  Blitz_LIBRARY_DIRS       = Location of Blitz libraries
#  Blitz_LIBRARIES          = Required libraries

######################################################################
#
# FindBlitz: find includes and libraries for BLITZ
#
# $Id: FindSciBlitz.cmake 1161 2011-12-17 15:44:00Z cary $
#
######################################################################

SciFindPackage(PACKAGE "Blitz"
              INSTALL_DIR "blitz"
              HEADERS "blitz/array.h"
              LIBRARIES "blitz"
              INCLUDE_SUBDIRS "include"
              LIBRARY_SUBDIRS "lib"
             )

if (BLITZ_FOUND)
  message(STATUS "Found Blitz")
  set(HAVE_BLITZ 1 CACHE BOOL "Whether have Blitz")
else ()
  message(STATUS "Did not find Blitz.  Use -DBLITZ_DIR to specify the installation directory.")
  if (SciBlitz_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
