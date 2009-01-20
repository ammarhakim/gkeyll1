##
# Find szip libs
##

from buildconf import BuildConf

bc = BuildConf(
# package name
    'SZIP',
# base paths for incs and libs
    ['$HOME/software/szip',
     '/usr/local/szip',
     '/usr/local',
     '/usr',
     ],
# headers to search for
    [],
# list of include directories
    [],
# libraries to look for
    ['szlib'],
# list of link directories
    ['lib']
    )
