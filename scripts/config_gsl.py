##
# Find GNU Scientific library (GSL)
##

from buildconf import BuildConf

bc = BuildConf(
# package name
    'GSL',
# base paths for incs and libs
    ['$HOME/software/gsl', '/usr'],
# headers to search for
    ['gsl/gsl_version.h', 'gsl/gsl_blas.h'],
# list of include directories
    ['include'],
# libraries to look for
    ['gsl', 'gslcblas'],
# list of link directories
    ['lib']
    )
