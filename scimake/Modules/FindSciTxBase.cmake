# - FindTxBase: Module to find include directories and
#   libraries for TxBase.
#
# Module usage:
#   find_package(TxBase ...)
#
# This module will define the following variables:
#  HAVE_TXBASE, TXBASE_FOUND = Whether libraries and includes are found
#  TxBase_INCLUDE_DIRS       = Location of TxBase includes
#  TxBase_LIBRARY_DIRS       = Location of TxBase libraries
#  TxBase_LIBRARIES          = Required libraries

######################################################################
#
# FindTxBase: find includes and libraries for txbase
#
# $Id: FindTxBase.cmake 1327 2012-04-23 17:45:30Z cary $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

if (NOT_HAVE_STD_ABS_DOUBLE)
  set(txbasefindlibs "txbase;txstd")
else ()
  set(txbasefindlibs txbase)
endif ()

if (ENABLE_PARALLEL)
  if (BUILD_WITH_SHARED_RUNTIME OR USE_SHARED_HDF5)
    set(instdirs txbase-parsh)
  else ()
    set(instdirs "txbase-par;txbase-ben")
  endif ()
else ()
  if (BUILD_WITH_SHARED_RUNTIME OR USE_SHARED_HDF5)
    set(instdirs txbase-sersh)
  else ()
    set(instdirs txbase)
  endif ()
endif ()

SciFindPackage(PACKAGE "TxBase"
  INSTALL_DIR "${instdirs}"
  HEADERS "txbase_version.h"
  LIBRARIES "${txbasefindlibs}"
  LIBRARY_SUBDIRS "lib/${CXX_COMP_LIB_SUBDIR};lib"
)

