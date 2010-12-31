######################################################################
#
# TxFortranChecks: check various C++ capabilities
#
# $Id: TxFortranChecks.cmake 167 2010-11-27 16:48:33Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

# Enable Fortran to all those variables
ENABLE_LANGUAGE(Fortran)
INCLUDE(FortranCInterface)

# Find the implicit link libraries
IF(DEBUG_CMAKE)
  MESSAGE("CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES = ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}")
  MESSAGE("CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES = ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
ENDIF()

# Remove mpi and system libs
SET(Fortran_IMPLICIT_LIBRARY_DIRS "")
SET(Fortran_IMPLICIT_LIBRARY_NAMES "")
SET(Fortran_IMPLICIT_LIBRARIES "")
SET(Fortran_IGNORED_LIBRARIES "")
FOREACH(txlib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  IF(${txlib} MATCHES ".*mpi.*")
    # MESSAGE("${txlib} is an MPI library.  Ignoring.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSEIF(${txlib} STREQUAL "opa" OR ${txlib} STREQUAL "dcmf.cnk" OR ${txlib} STREQUAL "dcmfcoll.cnk" OR ${txlib} STREQUAL "SPI.cna")
    # MESSAGE("${txlib} is a BGP compute node library.  Ignoring.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSEIF(${txlib} MATCHES "open-.*")
    # MESSAGE("${txlib} is an OpenMPI library.  Ignoring.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSEIF(${txlib} STREQUAL "dl" OR ${txlib} STREQUAL "nsl" OR ${txlib} STREQUAL "util" OR ${txlib} STREQUAL "m" OR ${txlib} STREQUAL "pthread" OR ${txlib} STREQUAL "c" OR ${txlib} STREQUAL "z")
    # MESSAGE("${txlib} is an system library.  Ignoring.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSEIF(${txlib} STREQUAL "sci_quadcore" OR ${txlib} STREQUAL "sma" OR ${txlib} STREQUAL "pmi" OR ${txlib} STREQUAL "rt" OR ${txlib} STREQUAL "alpslli" OR ${txlib} STREQUAL "alpsutil" OR ${txlib} STREQUAL "portals" OR ${txlib} STREQUAL "mv" OR ${txlib} STREQUAL "mpath" OR ${txlib} STREQUAL "pscrt")
    # MESSAGE("${txlib} is a cray pathscale library.  Not needed.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSEIF(${txlib} STREQUAL "pathfstart" OR ${txlib} STREQUAL "pathfortran" OR ${txlib} STREQUAL "z")
    # MESSAGE("${txlib} is a cray pathscale fortran library.  Not needed on CRAY.  May need to enable on some machines.")
    SET(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${txlib})
  ELSE()
    # MESSAGE("Looking for ${txlib}.")
    SET(txlibpathvar ${txlib}_LIBRARY)	# Cache variable
    FIND_LIBRARY(${txlibpathvar} ${txlib}
      PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
      NO_DEFAULT_PATH)
    SET(txlibpath ${${txlibpathvar}})
    IF(txlibpath)
      # MESSAGE("Adding ${txlibpath} to Fortran_IMPLICIT_LIBRARIES.")
      SET(Fortran_IMPLICIT_LIBRARIES ${Fortran_IMPLICIT_LIBRARIES} ${txlibpath})
      # MESSAGE("Adding ${txlib} to Fortran_IMPLICIT_LIBRARY_NAMES.")
      SET(Fortran_IMPLICIT_LIBRARY_NAMES ${Fortran_IMPLICIT_LIBRARY_NAMES} ${txlib})
      GET_FILENAME_COMPONENT(txlibdir ${txlibpath} PATH)
      # MESSAGE("Adding ${txlibdir} to Fortran_IMPLICIT_LIBRARY_DIRS.")
      SET(Fortran_IMPLICIT_LIBRARY_DIRS ${Fortran_IMPLICIT_LIBRARY_DIRS} ${txlibdir})
    ENDIF()
  ENDIF()
ENDFOREACH()
LIST(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARIES)
LIST(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARY_NAMES)
LIST(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARY_DIRS)
UNSET(Fortran_IMPLICIT_LIBFLAGS)
FOREACH(txlibdir ${Fortran_IMPLICIT_LIBRARY_DIRS})
  SET(Fortran_IMPLICIT_LIBFLAGS
    "${Fortran_IMPLICIT_LIBFLAGS} -L${txlibdir}")
ENDFOREACH()
FOREACH(txlib ${Fortran_IMPLICIT_LIBRARY_NAMES})
  SET(Fortran_IMPLICIT_LIBFLAGS
    "${Fortran_IMPLICIT_LIBFLAGS} -l${txlib}")
ENDFOREACH()
IF(Fortran_IMPLICIT_LIBFLAGS)
  STRING(STRIP ${Fortran_IMPLICIT_LIBFLAGS} Fortran_IMPLICIT_LIBFLAGS)
ENDIF()
IF(DEBUG_CMAKE)
  MESSAGE("Duplicates removed.  Result:")
  MESSAGE("Fortran_IMPLICIT_LIBRARIES = ${Fortran_IMPLICIT_LIBRARIES}")
  MESSAGE("Fortran_IMPLICIT_LIBRARY_NAMES = ${Fortran_IMPLICIT_LIBRARY_NAMES}")
  MESSAGE("Fortran_IMPLICIT_LIBRARY_DIRS = ${Fortran_IMPLICIT_LIBRARY_DIRS}")
  MESSAGE("Fortran_IMPLICIT_LIBFLAGS = ${Fortran_IMPLICIT_LIBFLAGS}")
  MESSAGE("Fortran_IGNORED_LIBRARIES = ${Fortran_IGNORED_LIBRARIES}")
ENDIF()

