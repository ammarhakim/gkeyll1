##
# Find MPI library
##

from wxbuildconf import WxBuildConf

bc = WxBuildConf(
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
