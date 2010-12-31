######################################################################
#
# TxFindCLapackCMake: find includes and libraries for txbase
#
# $Id: FindCLapackCMake.cmake 148 2010-10-22 10:15:38Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

TxFindPackage("CLapackCMake" "clapack_cmake"
  "clapack.h;f2c.h;blaswrap.h" "lapack;blas;f2c" "" "" "")

