##
# Find LAPACK library
##

import os
import re
from wxbuildconf import WxBuildConf

def makeLapackConfig(lapackLib):
    return WxBuildConf(
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
