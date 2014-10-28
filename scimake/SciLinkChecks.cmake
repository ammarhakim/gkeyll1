######################################################################
#
# SciLinkChecks: check/set various link flags
#
# $Id: SciLinkChecks.cmake 499 2014-02-08 23:41:38Z jrobcary $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

if (BUILD_SHARED_LIBS)
  if (WIN32)
# On windows define so that decl_spec stuff is done
    add_definitions(-DTXBASE_DLL)
  elseif (APPLE)
    if (NOT DEFINED CMAKE_INSTALL_NAME_DIR)
# Set library directory as needed by package managers
      set(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib)
    endif ()
  elseif (LINUX)
    set(CMAKE_EXE_LINKER_FLAGS "-Wl,-rpath,\$ORIGIN ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,-rpath,\$ORIGIN ${CMAKE_SHARED_LINKER_FLAGS}")
  endif ()
  if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
# Add the automatically determined parts of the RPATH that
# point to directories outside the build tree to the install RPATH
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif ()
  message(STATUS "After further modification:")
  foreach (type EXE SHARED)
    SciPrintVar(CMAKE_${type}_LINKER_FLAGS)
  endforeach ()
endif ()

