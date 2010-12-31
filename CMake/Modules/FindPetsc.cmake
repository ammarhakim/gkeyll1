######################################################################
#
# TxFindPetsc: find includes and libraries for petsc.  Complex
# due to the many libraries that petsc builds
#
# $Id: FindHdf5.cmake 148 2010-10-22 10:15:38Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

IF(ENABLE_PARALLEL)
  TxFindPackage("Petsc" "petsc-par" "petsc.h" 
      "petscts;petscsnes;petscksp;petscdm;petscmat;petscvec;petsc;cmumps;dmumps;smumps;zmumps;mumps_common;pord;scalapack;blacs;superlu_dist_2.3;HYPRE;parmetis;metis" 
      "" "" "")
ELSE()
  TxFindPackage("Petsc" "petsc" "petsc.h" 
      "petscts;petscsnes;petscksp;petscdm;petscmat;petscvec;petsc;cmumps;dmumps;smumps;zmumps;mumps_common;pord;scalapack;blacs;superlu_dist_2.3;HYPRE;parmetis;metis" 
      "" "" "")
ENDIF()

IF(PETSC_FOUND)
  SET(HAVE_PETSC 1 CACHE BOOL "Whether have the PETSC library")
ENDIF()

# The below should come from the petsc version.

