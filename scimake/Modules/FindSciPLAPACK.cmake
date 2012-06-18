# - FindSciGpulib: Module to find include directories and
#   libraries for Gpulib.
#
# Module usage:
#   find_package(SciGpulib ...)
#
# This module will define the following variables:
#  HAVE_GPULIB, GPULIB_FOUND = Whether libraries and includes are found
#  Gpulib_INCLUDE_DIRS       = Location of Gpulib includes
#  Gpulib_LIBRARY_DIRS       = Location of Gpulib libraries
#  Gpulib_LIBRARIES          = Required libraries

######################################################################
#
# FindGpulib.cmake: Find the GPU libraries
#
# $Id: FindSciPLAPACK.cmake 1342 2012-05-09 19:30:36Z lagae $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

SciFindPackage(
           PACKAGE "SciPLAPACK"
           INSTALL_DIR "txplapack"   # Empty?  Seems like we want a value here
           HEADERS PLA.h
           INCLUDE_SUBDIRS include
           LIBRARIES "txplapack"
           LIBRARY_SUBDIRS "lib/${CXX_COMP_LIB_SUBDIR};lib"
           )

if (SCIPLAPACK_FOUND)
  message(STATUS "Found txplapack")
  set(HAVE_SCIPLAPACK 1 CACHE BOOL "Whether have txplapack")
else ()
  message(STATUS "Did not find txplapack.  Use -DSCIPLAPACK_DIR to specify the installation directory.")
  if (SciPLAPACK_FIND_REQUIRED)
    message(FATAL_ERROR "Failed")
  endif ()
endif ()
