######################################################################
#
# Polyswift.cmake: Compute Polyswift specific options
#
# $Id: SciMpiLauncher.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

# First look in specified path
find_program(${sciexecvar}
  "${sciexec}"
  PATHS ${scipath}
  PATH_SUFFIXES "${sciexecsubdirs}"
  NO_DEFAULT_PATH
  DOC "Path to the ${sciexec} executable"
  )

# MPILAUNCHER for parallel runs
if (NOT DEFINED MPILAUNCHER)
    set(MPILAUNCHER ${MPILAUNCHER:FILEPATH})
endif ()

if (NOT DEFINED NPROCS)
  if (ENABLE_PARALLEL)
    set(NPROCS "2")
  endif ()
endif ()

