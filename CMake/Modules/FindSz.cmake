######################################################################
#
# TxFindSZ: find includes and libraries for sz (compression)
#
# $Id: FindSz.cmake 153 2010-10-27 22:59:15Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

TxFindPackage("SZ" "" "szlib.h" "sz" "sz" "" "")

# MESSAGE("Running FindSZ")

IF(SZ_FOUND)
  SET(HAVE_SZ 1 CACHE BOOL "Whether have the sz (compression) library")
ENDIF()

