##
# Find LAPACK library
##

import os
import re
from buildconf import BuildConf

def makeLapackConfig(lapackLib):
    return BuildConf(
        # package name
        'LAPACK',
        # base paths for incs and libs    
        [],
        # headers to search for
        [],
        # list of include directories
        [],
        # libraries to look for
        [lapackLib],
        # list of link directories
        ['$HOME/software/lapack',
         '/usr/lib/']
        )
