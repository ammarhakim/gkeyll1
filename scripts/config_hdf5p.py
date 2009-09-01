##
# Find HDF5 parallel library
##

import os
import re
from buildconf import BuildConf

bc = BuildConf(
# package name
    'HDF5P',
# base paths for incs and libs    
    ['$HOME/software/hdf5mpi',
     '/usr/local/hdf5mpi',
     '/usr/local',
     '/usr',
     ],
# headers to search for
    ['hdf5.h', 'H5public.h'],
# list of include directories
    ['include'],
# libraries to look for
    ['hdf5'],
# list of link directories
    ['lib']
    )

def testNewIfc(h5path):
    # determine the version of HDF5
    data = open(os.path.join(h5path, 'H5public.h')).read()
    majp = re.compile('#define\s+H5_VERS_MAJOR\s+(?P<num>\d)')
    minp = re.compile('#define\s+H5_VERS_MINOR\s+(?P<num>\d)')
    relp = re.compile('#define\s+H5_VERS_RELEASE\s+(?P<num>\d)')

    mo = majp.search(data)
    majv = int(mo.group('num'))
    mo = minp.search(data)
    minv = int(mo.group('num'))
    mo = relp.search(data)
    relv = int(mo.group('num'))

    new_ifc = False
    if majv > 1:
        new_ifc = True
    elif majv == 1:
        if minv > 6:
            new_ifc = True
        elif minv == 6:
            if relv >= 4:
                new_ifc = True

    return new_ifc
