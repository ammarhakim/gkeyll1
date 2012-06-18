# - FindSciTxGeom: Module to find include directories and
#   libraries for Geom.
#
# Module usage:
#   find_package(SciTxGeom ...)
#
# This module will define the following variables:
#  HAVE_TXGEOM, TXGEOM_FOUND = Whether libraries and includes are found
#  TxGeom_INCLUDE_DIRS       = Location of Geom includes
#  TxGeom_LIBRARY_DIRS       = Location of Geom libraries
#  TxGeom_LIBRARIES          = Required libraries


######################################################################
#
# FindSciTxGeom: find includes and libraries for txgeom
#
# $Id: FindSciTxGeom.cmake 1245 2012-01-31 21:36:22Z dws $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

if (HAVE_MPI)
  SciFindPackage(PACKAGE "TxGeom"
    INSTALL_DIR "txgeom-par;txgeom-ben"
    HEADERS "txgeom_version.h"
    LIBRARIES "txgeomtransforms;txgeomlibrary;txgeommesh;txgeomgeometry;txgeomtopology;txgeombase"
    LIBRARY_SUBDIRS "lib/${CXX_COMP_LIB_SUBDIR};lib"
  )
else ()
  SciFindPackage(PACKAGE "TxGeom"
    HEADERS "txgeom_version.h"
    LIBRARIES "txgeomtransforms;txgeomlibrary;txgeommesh;txgeomgeometry;txgeomtopology;txgeombase"
    LIBRARY_SUBDIRS "lib/${CXX_COMP_LIB_SUBDIR};lib"
  )
endif ()

if (TXGEOM_FOUND)
  # message(STATUS "Found TxGeom.")
  set(HAVE_TXGEOM 1 CACHE BOOL "Whether have TxGeom library")
else ()
  message(STATUS "Did not find TxGeom.  Use -DTXGEOM_DIR to specify the installation directory.")
  if (TxGeom_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()

