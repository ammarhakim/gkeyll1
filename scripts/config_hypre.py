##
# Find HYPRE library
##

import os
import re
from wxbuildconf import WxBuildConf

def makeHypreConfig(petscIncPath, petscLibPath, petscarch=None):
    r"""Returns an appropriate configure object for hypre library
    """
    if petscarch:
        incPath = '%s/include' % petscarch
        libPath = '%s/lib' % petscarch
    else:
        incPath = 'include'
        libPath = 'lib'
    
    return WxBuildConf(
        # package name
        'HYPRE',
        # base paths for incs and libs    
        ['$HOME/software/hypre',
         '/contrib/hypre',
         ],
        # headers to search for
        ['HYPRE.h'],
        # list of include directories
        [petscIncPath, incPath],
        # libraries to look for
        ['HYPRE'],
        # list of link directories
        [petscLibPath, libPath]
        )
