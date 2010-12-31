######################################################################
#
# TxFindTrilinos: find includes and libraries for Trilinos
#
# $Id: FindTrilinos.cmake 148 2010-10-22 10:15:38Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

IF(ENABLE_PARALLEL)
  IF(WIN32)
    TxFindPackage("Trilinos" "trilinos-par" "az_aztec.h"
      "aztecoo;ml;amesos;ifpack;epetraext;galeri;triutils;epetra;teuchos"
      "" "" "")
  ELSE()
    TxFindPackage("Trilinos" "trilinos-par" "az_aztec.h"
      "aztecoo;ml;zoltan;amesos;ifpack;epetraext;galeri;triutils;epetra;teuchos"
      "" "" "")
  ENDIF()
ELSE()
  TxFindPackage("Trilinos" "trilinos" "az_aztec.h"
      "aztecoo;ml;amesos;ifpack;epetraext;galeri;triutils;epetra;teuchos"
      "" "" "")
ENDIF()

IF(TRILINOS_FOUND)
  OPTION(HAVE_AMESOS "Trilinos Amesos Library" ON)
  OPTION(HAVE_EPETRAEXT "Trilinos Epetraext library" ON)
  OPTION(HAVE_TRILINOS "Trilinos libraries" ON)

# Find the blas and lapack used by Trilinos.
  FIND_FILE(TRILINOS_CONFIG_CMAKE
    TrilinosConfig.cmake
    PATHS ${Trilinos_INCLUDE_DIR} ${Trilinos_INCLUDE_DIR}/..
  )
  IF(NOT TRILINOS_CONFIG_CMAKE)
    MESSAGE(FATAL_ERROR "TRILINOS_CONFIG_CMAKE not found.")
  ENDIF()
  IF(FALSE)
    INCLUDE(${TRILINOS_CONFIG_CMAKE})
  ELSE()
    FILE(STRINGS ${TRILINOS_CONFIG_CMAKE} tritpls
      REGEX ".*Trilinos_TPL_LIBRARIES.*")
    # MESSAGE("tritpls = ${tritpls}")
    STRING(REGEX REPLACE "SET.Trilinos_TPL_LIBRARIES ." "" tritpls "${tritpls}")
    # MESSAGE("tritpls = ${tritpls}")
    STRING(REGEX REPLACE "..$" "" tritpls "${tritpls}")
    # MESSAGE("tritpls = ${tritpls}")
    STRING(REGEX REPLACE ".;" ";" Trilinos_TPL_LIBRARIES "${tritpls}")
# Make this a string list
    SEPARATE_ARGUMENTS(Trilinos_TPL_LIBRARIES)
# with no duplicates
    LIST(REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
    MESSAGE("Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  ENDIF()
ENDIF(TRILINOS_FOUND)

