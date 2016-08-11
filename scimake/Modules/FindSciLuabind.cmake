# - FindSciLuabind: Module to find include directories and
#   libraries for Luabind.
#
# Module usage:
#   find_package(SciLuabind ...)
#
# This module will define the following variables:
#  HAVE_LUABIND, LUABIND_FOUND = Whether libraries and includes are found
#  Luabind_INCLUDE_DIRS       = Location of Luabind includes
#  Luabind_LIBRARY_DIRS       = Location of Luabind libraries
#  Luabind_LIBRARIES          = Required libraries

######################################################################
#
# FindLuabind: find includes and libraries for LUABIND
#
# $Id: FindSciLuabind.cmake 1161 2011-12-17 15:44:00Z cary $
#
######################################################################

SciFindPackage(PACKAGE "Luabind"
              INSTALL_DIR "luabind"
              HEADERS "luabind/luabind.hpp"
              LIBRARIES "luabind"
              INCLUDE_SUBDIRS "include"
              LIBRARY_SUBDIRS "lib"
             )

if (LUABIND_FOUND)
  message(STATUS "Found Luabind")
  set(HAVE_LUABIND 1 CACHE BOOL "Whether have Luabind")
else ()
  message(STATUS "Did not find Luabind.  Use -DLUABIND_DIR to specify the installation directory.")
  if (SciLuabind_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
