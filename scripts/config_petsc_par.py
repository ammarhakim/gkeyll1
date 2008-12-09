##
# Find PETSC parallel library
##

import os
import re
from buildconf import BuildConf

def makePetscConfig(petscarch):
    r"""Returns an appropriate configure object for PETSC PARALLEL
    """
    if petscarch:
        return BuildConf(
            # package name
            'PETSCPAR',
            # base paths for incs and libs
            ['$HOME/software/petscmpi',
             '/contrib/petscmpi',
             '/usr/local/petscmpi',
             ],
            # headers to search for
            ['petsc.h', 'petscconf.h'],
            # list of include directories
            ['include'],
            # libraries to look for
            ['petsccontrib', 'petscts', 'petscsnes', 'petscksp', 'petscdm', \
             'petscmat', 'petscvec', 'petsc'],
            # list of link directories
            ['lib/%s' % petscarch],
            # list of sub-includes to add to include paths
            ['../bmake/%s' % petscarch]
            )
    else:
        return BuildConf(
            # package name
            'PETSCPAR',
            # base paths for incs and libs
            ['$HOME/software/petscmpi',
             '/contrib/petscmpi',
             '/usr/local/petscmpi',
             ],
            # headers to search for
            ['petsc.h', 'petscconf.h'],
            # list of include directories
            ['include'],
            # libraries to look for
            ['petsccontrib', 'petscts', 'petscsnes', 'petscksp', 'petscdm', \
             'petscmat', 'petscvec', 'petsc'],
            # list of link directories
            ['lib']
            )
