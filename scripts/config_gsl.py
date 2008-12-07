##
# Find GNU Scientific library (GSL)
##

from wxbuildconf import WxBuildConf

bc = WxBuildConf(
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
