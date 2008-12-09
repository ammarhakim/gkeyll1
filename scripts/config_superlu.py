##
# Find SUPERLU library
##

import os
import re
from buildconf import BuildConf

def makeSuperLUConfig(petscIncPath, petscLibPath, petscarch=None):
    r"""Returns an appropriate configure object for superlu library
    """
    
    if petscarch:
        incPath = '%s/SRC' % petscarch
        libPath = '%s/' % petscarch
    else:
        incPath = 'include'
        libPath = 'lib'
    
    return BuildConf(
        # package name
        'SUPERLU',
        # base paths for incs and libs    
        ['$HOME/software/superlu',
         '/contrib/superlu',
         ],
        # headers to search for
        ['colamd.h', 'old_colamd.h' ,'slu_Cnames.h', 'slu_ddefs.h'], #, 'slu_util.h'],
        # list of include directories
        [petscIncPath, incPath],
        # libraries to look for
        ['superlu_3.0'],
        # list of link directories
        [petscLibPath, libPath]
        )
