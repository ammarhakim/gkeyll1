# - FindSciMagma: Module to find the MAGMA library.
#
# Module usage:
#   find_package(SciMagma ...)
#
# Should probably be modified to use SciFindPackage...

######################################################################
#
# FindSciMagma.cmake: Find the MAGMA library.
#
# $Id: FindSciMagma.cmake 342 2013-08-21 17:17:10Z michaelgalloy $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

if (WIN32)
  set(MAGMA_LIB_PREFIX "lib")
  set(MAGMA_LIB_SUFFIX ".lib")
else ()
  set(MAGMA_LIB_PREFIX "lib")
  set(MAGMA_LIB_SUFFIX ".a")
endif ()

SciFindPackage(PACKAGE "magma"
  INSTALL_DIR "${magma_ROOT_DIR}"
  LIBRARIES "${MAGMA_LIB_PREFIX}magma${MAGMA_LIB_SUFFIX}"
  HEADERS "magma.h"
  INCLUDE_SUBDIRS "include"
  LIBRARY_SUBDIRS "lib"
  ALLOW_LIBRARY_DUPLICATES TRUE
)

if (MAGMA_FOUND)
  message(STATUS "Found MAGMA.")
  set(HAVE_MAGMA 1 CACHE BOOL "Whether have MAGMA")
else ()
   if (SciMagma_FIND_REQUIRED)
     message(FATAL_ERROR "Could not find MAGMA")
   else ()
     if (magma_ROOT_DIR)
       message(STATUS "MAGMA not found in ${magma_ROOT_DIR}")
     else ()
       message(STATUS "Not searching for MAGMA")
     endif ()
   endif ()
endif ()
