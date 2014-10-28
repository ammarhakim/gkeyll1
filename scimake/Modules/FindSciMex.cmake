# - FindSciMex: Module to find include directories and
#   libraries for Mex.
#
# Module usage:
#   find_package(SciMex ...)
#
# May need to be changed to use SciFindPackge()

########################################################################
#
# FindSciMex
#
# $Id: FindSciMex.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
########################################################################

find_file(MEX NAMES mex.bat mex HINTS ENV PATH)

if (MEX)
  set(MEX_FOUND TRUE)
endif ()

if (MEX_FOUND)
  if (NOT MEX_FIND_QUIETLY)
    message(STATUS "Found MEX: ${MEX}")
  endif ()
  set(HAVE_MEX 1 CACHE BOOL "Whether have MEX")
else ()
   if (SciMex_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find MEX")
   endif ()
endif ()
