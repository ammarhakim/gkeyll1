##
# Find PETSC serial library
##

import os
import re
from wxbuildconf import WxBuildConf

def makePetscConfig(petscarch):
    r"""Returns an appropriate configure object for PETSC
    """
    if petscarch:
        return WxBuildConf(
            # package name
            'PETSCSER',
            # base paths for incs and libs
            ['$HOME/software/petsc',
             '/contrib/petsc',
             '/usr/local/petsc',
             ],
            # headers to search for
            ['petsc.h', 'petscconf.h'],
            # list of include directories
            ['include'],
            # libraries to look for
            ['petsccontrib', 'petscts', 'petscsnes', 'petscksp', 'petscdm', \
             'petscmat', 'petscvec', 'petsc', 'mpiuni'],
            # list of link directories
            ['lib/%s' % petscarch],
            # list of sub-includes to add to include paths
            ['mpiuni', '../bmake/%s' % petscarch]
            )
    else:
        return WxBuildConf(
            # package name
            'PETSCSER',
            # base paths for incs and libs
            ['$HOME/software/petsc',
             '/contrib/petsc',
             '/usr/local/petsc',
             ],
            # headers to search for
            ['petsc.h', 'petscconf.h'],
            # list of include directories
            ['include'],
            # libraries to look for
            ['petsccontrib', 'petscts', 'petscsnes', 'petscksp', 'petscdm', \
             'petscmat', 'petscvec', 'petsc', 'mpiuni'],
            # list of link directories
            ['lib'],
            # list of sub-includes to add to include paths
            ['mpiuni']
            )
