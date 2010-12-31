######################################################################
#
# FindTxBlasLapack: find includes and libraries for lapack and blas
#
# $Id: FindTxBlasLapack.cmake 87 2010-07-23 18:08:52Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# The goal is to find BLAS and LAPACK, one implementation
# or another.  It should be invoked after Trilinos or PETSC
# have been found, if they are used.
#
# User specified invocation args, like
#  -DLAPACK_LIBRARY_DIRS:PATH='<lapackdir>' -DLAPACK_LIBRARY_NAMES:STRING='lapack' -DBLAS_LIBRARY_DIRS:PATH='<blasdir>' -DBLAS_LIBRARY_NAMES:STRING='blas'
# are determining.  If not there, a search is done.
#
# Options allow one to first search by
#   using whatever trilinos linked against
#   using whatever PETSC linked against
#
# User specified.
# Optionally: use whatever trilinos linked against
#   if trilinos found.  (This is the VORPAL preference.)
# Optionally: use whatever PETSc linked against if PETSc
#   found.  (This is the FACETS preference.)
# System optimized (use FCLIBS or equiv to help link):
#   -framework Accelerate
#   essl
#   mkl
#   acml
#   /contrib/atlas-ser (the one built againt lapack)
#   /contrib/altas-clp (the one built against clapack)
#   Other ATLAS in system directories
#   Other system lapack-blas in system directories.
#   Optionally GOTO (as it has a funky license)
# lapack in /contrib
# clapack in /contrib
#
######################################################################

# TxFindPackage("CLapackCMake" "clapack_cmake"
#  "clapack.h;f2c.h;blaswrap.h" "lapack;blas;f2c" ""
#  "" "" "" "NO_DEFAULT_PATH")

