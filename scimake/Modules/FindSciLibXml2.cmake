# - FindSciLibXml2n: Module to find include directories and libraries
#   for LibXml2. This module was implemented as there is no stock
#   CMake module for LibXml2. 
#
# This module can be included in CMake builds in find_package:
#   find_package(SciLibXml2 REQUIRED)
#
# This module will define the following variables:
#  HAVE_LIBXML2         = Whether have the LibXml2 library
#  LibXml2_INCLUDE_DIRS = Location of LibXml2 includes
#  LibXml2_LIBRARY_DIRS = Location of LibXml2 libraries
#  LibXml2_LIBRARIES    = Required libraries
#  LibXml2_STLIBS       = Location of LibXml2 static library

######################################################################
#
# SciFindLibXml2: find includes and libraries for LibXml2.
#
# $Id: FindSciLibXml2.cmake 1323 2012-04-19 18:59:43Z dws $
#
# Copyright 2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################
set(SUPRA_SEARCH_PATH ${SUPRA_SEARCH_PATH})

SciFindPackage(PACKAGE "LibXml2"
              INSTALL_DIR "libxml2"
              HEADERS "include/libxml2"
              LIBRARIES "libxml2.a"
              )

if (LIBXML2_FOUND)
  message(STATUS "Found LibXml2")
  set(HAVE_LIBXML2 1 CACHE BOOL "Whether have the LIBXML2 library")
else ()
  message(STATUS "Did not find LibXml2.  Use -DLIBXML2_DIR to specify the installation directory.")
  if (SciLibXml2_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()

