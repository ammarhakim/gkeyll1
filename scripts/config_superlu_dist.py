##
# Find SUPERLU library
##

import os
import re
from wxbuildconf import WxBuildConf

def makeSuperLUDistConfig(petscIncPath, petscLibPath, petscarch=None):
    r"""Returns an appropriate configure object for superlu library
    """
    
    if petscarch:
        incPath = '%s/SRC' % petscarch
        libPath = '%s/' % petscarch
    else:
        incPath = 'include'
        libPath = 'lib'
    
    return WxBuildConf(
        # package name
        'SUPERLUDIST',
        # base paths for incs and libs    
        ['$HOME/software/superlu',
         '/contrib/superlu',
         ],
        # headers to search for
        ['superlu_defs.h','superlu_ddefs.h','superlu_zdefs.h',\
        'util_dist.h','machines.h','Cnames.h','dcomplex.h','old_colamd.h'], #, 'slu_util.h'],
        # list of include directories
        [petscIncPath, incPath],
        # libraries to look for
        ['superlu_dist_2.0'],
        # list of link directories
        [petscLibPath, libPath]
        )
