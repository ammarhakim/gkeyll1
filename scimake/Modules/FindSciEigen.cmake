# - FindSciEigen: Module to find include directories and
#   libraries for Eigen.
#
# Module usage:
#   find_package(SciEigen ...)
#
# This module will define the following variables:
#  HAVE_EIGEN, EIGEN_FOUND = Whether libraries and includes are found
#  Eigen_INCLUDE_DIRS       = Location of Eigen includes
#  Eigen_LIBRARY_DIRS       = Location of Eigen libraries
#  Eigen_LIBRARIES          = Required libraries

######################################################################
#
# FindEigen: find includes and libraries for EIGEN
#
# $Id: FindSciEigen.cmake 1161 2011-12-17 15:44:00Z cary $
#
######################################################################

SciFindPackage(PACKAGE "Eigen"
              INSTALL_DIR "eigen"
              HEADERS "Eigen/Core"
              LIBRARIES ""
              INCLUDE_SUBDIRS "include/eigen3"
              LIBRARY_SUBDIRS ""
             )

if (EIGEN_FOUND)
  message(STATUS "Found Eigen")
  set(HAVE_EIGEN 1 CACHE BOOL "Whether have Eigen")
else ()
  message(STATUS "Did not find Eigen.  Use -DEIGEN_DIR to specify the installation directory.")
  if (SciEigen_FIND_REQUIRED)
    message(FATAL_ERROR "Failed.")
  endif ()
endif ()
