##
# Find BLAS library
##

import os
import re
from wxbuildconf import WxBuildConf

def makeBlasConfig(blasLib):
    return WxBuildConf(
        # package name
        'BLAS',
        # base paths for incs and libs    
        [],
        # headers to search for
        [],
        # list of include directories
        [],
        # libraries to look for
        [blasLib],
        # list of link directories
        ['$HOME/software/blas',
         '/usr/lib/']
        )
