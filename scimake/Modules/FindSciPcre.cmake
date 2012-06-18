######################################################################
#
# FindPCRE: find includes and libraries for PCRE
#
# $Id: FindSciPcre.cmake 1324 2012-04-21 13:31:24Z cary $
#
######################################################################

# - FindSciPCRE: Module to find include directories and
#   libraries for PCRE
#
# Module usage:
#   find_package(SciPcre ...)
#
# This module will define the following variables:
#  HAVE_PCRE, PCRE_FOUND  = Whether libraries and includes are found
#  Pcre_INCLUDE_DIRS       = Location of Gsl includes
#  Pcre_LIBRARY_DIRS       = Location of Gsl libraries
#  Pcre_LIBRARIES          = Required libraries

SciFindPackage(PACKAGE "Pcre"
              INSTALL_DIR "pcre"
              HEADERS "pcre.h"
              LIBRARIES "pcre"
              INCLUDE_SUBDIRS "include"
              LIBRARY_SUBDIRS "lib"
             )

if (PCRE_FOUND)
  message(STATUS "Found Pcre")
  set(HAVE_PCRE 1 CACHE BOOL "Whether have Pcre")
else ()
  message(STATUS "Did not find Pcre.  Use -DPCRE_DIR to specify the installation directory.")
  if (SciPcre_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
