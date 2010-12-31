######################################################################
#
# FindGsl: find includes and libraries for txgeom
#
# $Id: FindGsl.cmake 144 2010-10-19 07:01:22Z cary $
#
######################################################################

TxFindPackage("Gsl" "gsl" "gsl/gsl_math.h;gsl/gsl_sf_legendre.h"
  "gsl;gslcblas" "" "" "lib" "" "")
