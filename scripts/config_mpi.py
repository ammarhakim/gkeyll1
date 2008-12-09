##
# Find MPI library
##

from buildconf import BuildConf

bc = BuildConf(
# package name
    'MPI',
# base paths for incs and libs    
    ['/usr/local/mpi',
     '/usr',
     ],
# headers to search for
    [],
# list of include directories
    [],
# libraries to look for
    [],
# list of link directories
    []
    )
# also add mpicxx for binaries to search for
bc.setBin('mpicxx', ['bin'])
