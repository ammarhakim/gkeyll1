# - FindSciNsis: Module to find include directories and
#   libraries for Nsis.
#
# Module usage:
#   find_package(SciNsis ...)
#
# May need to be changed to use SciFindPackage()

#################################################
#
# Find NSIS packager
#
# $Id: FindSciNsis.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
#################################################

if (WIN32)
  set(SCIC_NSIS_SEARCHPATH
    "$ENV{PROGRAMFILES}/NSIS" "$ENV{PROGRAMFILES(X86)}/NSIS"
  )
endif ()

find_program(MAKENSIS
  makensis
  PATHS ${SCIC_NSIS_SEARCHPATH}
  DOC "Location of the NSIS executable"
)
if (MAKENSIS)
  set(MAKENSIS_FOUND TRUE)
endif ()

