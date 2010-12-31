######################################################################
#
# TxFindMpi: check whether the compiler wraps MPI, if not, find MPI
#
# $Id: TxFindMpi.cmake 179 2010-12-15 21:17:05Z dws $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

OPTION(ENABLE_PARALLEL "Enable parallel build")

IF(ENABLE_PARALLEL)
  MESSAGE(STATUS "ENABLE_PARALLEL requested.")
# See whether mpi is already in the compiler
  TRY_COMPILE(TX_HAVE_MPI_COMPILER_WRAPPER ${PROJECT_BINARY_DIR}/CMake
    ${PROJECT_SOURCE_DIR}/CMake/mpi_h.cpp)
  IF(TX_HAVE_MPI_COMPILER_WRAPPER)
    SET(MPI_FOUND TRUE)
    MESSAGE(STATUS "${PROJECT_SOURCE_DIR}/CMake/mpi_h.cpp compiled.")
  ELSE()
    MESSAGE(STATUS "${PROJECT_SOURCE_DIR}/CMake/mpi_h.cpp did not compile.")
    FIND_PACKAGE(MPI REQUIRED)
    IF(MPI_FOUND)
      INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
      MESSAGE(STATUS "mpi include path is: ${MPI_INCLUDE_PATH}")
    ELSE(MPI_FOUND)
      MESSAGE(FATAL_ERROR "Parallel requested but MPI not found.")
    ENDIF(MPI_FOUND)
  ENDIF()
ELSE(ENABLE_PARALLEL)
  MESSAGE(STATUS "ENABLE_PARALLEL not requested.")
ENDIF(ENABLE_PARALLEL)

IF(MPI_FOUND)
  SET(HAVE_MPI 1 CACHE BOOL "Whether using MPI")
  MESSAGE(STATUS "Enabling MPI")
#  MESSAGE(WARNING "Make sure library dependencies (such as hdf5 and Trilinos) are the MPI versions")
ELSE()
  MESSAGE(STATUS "MPI not enabled.")
ENDIF()

