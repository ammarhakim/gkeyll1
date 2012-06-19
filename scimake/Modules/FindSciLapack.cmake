# - FindSciCLapackscimake: Module to find include directories and
#   libraries for CLapackscimake.
#
# Module usage:
#   find_package(SciCLapackscimake ...)
#
# This module will define the following variables:
#  HAVE_LAPACK, LAPACK_FOUND = Whether libraries and includes are found
#  LAPACK_INCLUDE_DIRS = Location of CLapackscimake includes
#  LAPACK_LIBRARY_DIRS = Location of CLapackscimake libraries
#  LAPACK_LIBRARIES    = Required libraries

######################################################################
#
# FindLapack: find includes and libraries for txbase
#
# $Id: FindSciLapack.cmake 1245 2012-01-31 21:36:22Z dws $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

SciFindPackage(PACKAGE "Lapack"
              INSTALL_DIR "lapack"
              HEADERS ""
              LIBRARIES "lapack;blas"
              )

if (LAPACK_FOUND)
  message(STATUS "Lapack found.")
  set(HAVE_LAPACK 1 CACHE BOOL "Whether have Lapack")
else ()
  message(STATUS "Did not find Lapack.  Use -DLAPACKE_DIR to specify the installation directory.")
  if (SciCLapackscimake_FIND_REQUIRED)
    message(FATAL_ERROR "Failed")
  endif ()
endif ()
