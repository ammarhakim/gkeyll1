# - FindSciMkl: Module to find include directories and
#   libraries for Mkl.
#
# Module usage:
#   find_package(SciMkl ...)
#
# This module will define the following variables:
#  HAVE_MKL, MKL_FOUND = Whether libraries and includes are found
#  Mkl_INCLUDE_DIRS = Location of Mkl includes
#  Mkl_LIBRARY_DIRS = Location of Mkl libraries
#  Mkl_LIBRARIES    = Required libraries
#  Mkl_STLIB        = Static libraries
#  Iomp5_LIBRARIES  = Openmp intel libraries

######################################################################
#
# FindSciMkl: find includes and libraries for txbase
#
# $Id: FindSciMkl.cmake 484 2014-01-26 16:39:04Z jrobcary $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

if (WIN32)
  #foreach(year 2011 2012 2013)
    #set(Mkl_ROOT_DIR "C:/Program Files (x86)/Intel/Composer XE ${year}/mkl/lib/intel64")
    set(Mkl_ROOT_DIR "C:/Program Files (x86)/Intel/Composer XE/mkl")
    set(Mkl_INCLUDE_DIRS "${Mkl_ROOT_DIR}/include")
      SciFindPackage(PACKAGE "Mkl"
                    LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core"
                    INCLUDE_SUBDIRS "include"
                    LIBRARY_SUBDIRS "lib/intel64"
                    )
    #set(Iomp5_ROOT_DIR "C:/Program Files (x86)/Intel/Composer XE ${year}/compiler/lib/intel64")
    set(Iomp5_ROOT_DIR "C:/Program Files (x86)/Intel/Composer XE/compiler/lib/intel64")
      SciFindPackage(PACKAGE "Iomp5"
                    LIBRARIES "libiomp5md"
                    )
    if (MKL_FOUND)
      message(STATUS "Mkl found.")
      set(HAVE_MKL 1 CACHE BOOL "Whether have Mkl")
      break()
    endif ()
  #endforeach()
else (WIN32)
    if (NOT Mkl_ROOT_DIR)
      set(Mkl_ROOT_DIR "/usr/local/intel/mkl")
    endif ()
    if (NOT Iomp5_ROOT_DIR)
      set (Iomp5_ROOT_DIR "${Mkl_ROOT_DIR}/../lib/intel64")
    endif ()
    set(Mkl_INCLUDE_DIRS "${Mkl_ROOT_DIR}/include")
      if (NOT LINK_WITH_MKL_SHARED)
        SciFindPackage(PACKAGE "Mkl"
                       LIBRARIES "libmkl_intel_lp64.a;libmkl_intel_thread.a;libmkl_core.a"
                       INCLUDE_SUBDIRS "include"
                       LIBRARY_SUBDIRS "lib/intel64"
                       )
      else ()
         # On iter one needs to link against the shared MKL libraries
         SciFindPackage(PACKAGE "Mkl"
                       LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core"
                       INCLUDE_SUBDIRS "include"
                       LIBRARY_SUBDIRS "lib/intel64"
                       )

      endif ()
    SciFindPackage(PACKAGE "Iomp5"
                   LIBRARIES "iomp5"
                   )
    if (MKL_FOUND)
      message(STATUS "Mkl found.")
      set(HAVE_MKL 1 CACHE BOOL "Whether have Mkl")
    endif ()
endif (WIN32)

if (NOT MKL_FOUND)
  message(STATUS "Did not find Mkl.  Use -DMkl_ROOT_DIR and, if using OpenMP, Iomp5_ROOT_DIR to specify the installation directory.")
  if (SciMkl_FIND_REQUIRED)
    message(FATAL_ERROR "Finding MKL failed.")
  endif ()
endif ()

