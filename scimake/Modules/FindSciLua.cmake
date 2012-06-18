# - FindSciLua: Module to find include directories and
#   libraries for Lua.
#
# Module usage:
#   find_package(SciLua ...)
#
# This module will define the following variables:
#  HAVE_LUA, LUA_FOUND = Whether libraries and includes are found
#  Lua_INCLUDE_DIRS       = Location of Lua includes
#  Lua_LIBRARY_DIRS       = Location of Lua libraries
#  Lua_LIBRARIES          = Required libraries

######################################################################
#
# FindLua: find includes and libraries for LUA
#
# $Id: FindSciLua.cmake 1161 2011-12-17 15:44:00Z cary $
#
######################################################################

SciFindPackage(PACKAGE "Lua"
              INSTALL_DIR "lua"
              HEADERS "lauxlib.h;lua.h;lua.hpp;luaconf.h;lualib.h"
              LIBRARIES "lua"
              INCLUDE_SUBDIRS "include"
              LIBRARY_SUBDIRS "lib"
             )

if (LUA_FOUND)
  message(STATUS "Found Lua")
  set(HAVE_LUA 1 CACHE BOOL "Whether have Lua")
else ()
  message(STATUS "Did not find Lua.  Use -DLUA_DIR to specify the installation directory.")
  if (SciLua_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
