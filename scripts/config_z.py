##
# Find Z libs
##

import os
from buildconf import BuildConf

if os.name == 'nt':
    zlib = 'zlib'
else:
    zlib = 'z'

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
    [zlib],
# list of link directories
    ['lib']
    )
