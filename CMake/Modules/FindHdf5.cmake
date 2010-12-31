######################################################################
#
# TxFindHdf5: find includes and libraries for txbase
#
# $Id: FindHdf5.cmake 148 2010-10-22 10:15:38Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

IF(ENABLE_PARALLEL)
  TxFindPackage("Hdf5" "hdf5-par" "hdf5.h" "hdf5" "hdf5" "" "")
ELSE()
  TxFindPackage("Hdf5" "hdf5" "hdf5.h" "hdf5" "hdf5" "" "")
ENDIF()

IF(HDF5_FOUND)
  SET(HAVE_HDF5 1 CACHE BOOL "Whether have the HDF5 library")
  SET(OLD_H5S_SELECT_HYPERSLAB_IFC 0 CACHE BOOL
    "Whether using the old 1.6.3 H5Sselect_hyperslab interface")
ENDIF()

# The below should come from the hdf5 version.

