######################################################################
#
# TxFindZ: find includes and libraries for z (compression)
#
# $Id: FindZ.cmake 148 2010-10-22 10:15:38Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

TxFindPackage("Z" "" "zlib.h" "z" "z" "" "")

IF(Z_FOUND)
  SET(HAVE_Z 1 CACHE BOOL "Whether have the z (compression) library")
ENDIF()

