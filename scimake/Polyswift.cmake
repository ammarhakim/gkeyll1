######################################################################
#
# Polyswift.cmake: Compute Polyswift specific options
#
# $Id: Polyswift.cmake 1245 2012-01-31 21:36:22Z dws $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

# Set allowed number of spatial dimensions
#IF(NOT DEFINED POLYSWIFT_NDIMS)
#  SET(POLYSWIFT_NDIMS "1;2;3")
#ENDIF()

#MESSAGE(STATUS "Enabling NDIMS of ${POLYSWIFT_NDIMS}")

#SET(ndimmax 0)
#FOREACH(ndim ${POLYSWIFT_NDIMS})
#  SET(NDIM${ndim} 1 CACHE BOOL "Define to compile NDIM = ${ndim}")
#  IF(${ndim} GREATER ${ndimmax})
#    SET(ndimmax ${ndim})
#  ENDIF()
#ENDFOREACH()
#SET(NDIMMAX ${ndimmax} CACHE STRING "Maximum value of compiled NDIM")

