##
# Find Z libs
##

from buildconf import BuildConf

bc = BuildConf(
# package name
    'ZLIB',
# base paths for incs and libs
    ['$HOME/software/zlib',
     '/usr/local/zlib',
     '/usr/local',
     '/usr',
     ],
# headers to search for
    [],
# list of include directories
    [],
# libraries to look for
    ['z'],
# list of link directories
    ['lib']
    )
