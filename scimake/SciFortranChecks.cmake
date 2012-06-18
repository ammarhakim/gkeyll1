######################################################################
#
# SciFortranChecks: check various Fortran capabilities
#
# $Id: SciFortranChecks.cmake 1357 2012-06-05 19:50:29Z cary $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# THIS FILE NEEDS ALOT OF WORK.  JUST STARTING WITH GNU.
#
######################################################################

include(${SCICMAKE_DIR}/SciFortranFindVersion.cmake)

# Set the lib subdir from the Compiler ID and version
if (DEBUG_CMAKE)
  SciPrintString("CMAKE_Fortran_COMPILER_ID = ${CMAKE_Fortran_COMPILER_ID}.")
endif ()
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL GNU)
  # Not sure if we need an option to specify whether or not to add
  # openmp: Just always add it for now.  Another alternative is to add
  # an openmp variable that can be added on a case-by-case basis.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
  if (NOT USING_MINGW)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pipe")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fPIC -O3")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fPIC -O0 -g -Wp,-DDEBUG")
  endif ()
  string(SUBSTRING ${Fortran_VERSION} 0 3 Fortran_MAJOR_VERSION)
  set(Fortran_COMP_LIB_SUBDIR gfortran${Fortran_MAJOR_VERSION})
  set(FC_DOUBLE_FLAGS "-fdefault-real-8 -fdefault-double-8")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL Intel)
  string(REGEX REPLACE "\\.[0-9]+.*$" "" Fortran_MAJOR_VERSION ${Fortran_VERSION})
  set(Fortran_COMP_LIB_SUBDIR icpc${Fortran_MAJOR_VERSION})
  set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} iomp5)
  if (WIN32)
    set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -NODEFAULTLIB:MSVCRT")
  endif ()
  set(FC_DOUBLE_FLAGS "-autodouble")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL PathScale)
  string(SUBSTRING ${Fortran_VERSION} 0 1 Fortran_MAJOR_VERSION)
  set(Fortran_COMP_LIB_SUBDIR path${Fortran_MAJOR_VERSION})
  set(FC_DOUBLE_FLAGS "-r8")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL PGI)
  string(REGEX REPLACE "\\.[0-9]+-.*$" "" Fortran_MAJOR_VERSION ${Fortran_VERSION})
  set(Fortran_COMP_LIB_SUBDIR pgi${Fortran_MAJOR_VERSION})
# Compiler optimization flags set based on "ultra" optimization in
# flags.m4.  Overrides scimake default, since that had -Mipa=fast (no inline).
  set(CMAKE_Fortran_FLAGS_RELEASE
    "-fast -O3 -DNDEBUG -Munroll -Minline=levels:5 -Mipa=fast,inline -Mmovnt")
# For a fully-optimized build, set IPA options for linker too
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE
    "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Mipa=fast,inline")
  set(FC_DOUBLE_FLAGS "-r8")
elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL XL)
# This should be the basename of the compiler
  string(REGEX REPLACE "\\.[0-9]+.*$" "" Fortran_MAJOR_VERSION ${Fortran_VERSION})
  string(REGEX REPLACE "^0+" "" Fortran_MAJOR_VERSION ${Fortran_MAJOR_VERSION})
#SriV. want file name without extension.
  get_filename_component(REL_CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER} NAME_WE)
# Since we install ben builds in a completely different directory, can
# use same name for Fortran_COMP_LIB_SUBDIR
  if (${REL_CMAKE_Fortran_COMPILER} MATCHES ".*_r$")
    set(Fortran_COMP_LIB_SUBDIR xlC_r${Fortran_MAJOR_VERSION})
  else ()
    set(Fortran_COMP_LIB_SUBDIR xlC${Fortran_MAJOR_VERSION})
  endif ()
  set(SEPARATE_INSTANTIATIONS 1 CACHE BOOL "Whether to separate instantiations -- for correct compilation on xl")
  set(FC_DOUBLE_FLAGS "-qrealsize=8 -qautodbl=dbl4")
endif ()
if (SCI_FC_PROMOTE_REAL_TO_DOUBLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FC_DOUBLE_FLAGS}")
endif ()
SciPrintString("Fortran_COMP_LIB_SUBDIR = ${Fortran_COMP_LIB_SUBDIR}")

SciPrintString("RESULTS FOR cmake detected fortran implicit libraries:")
SciPrintVar(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
SciPrintVar(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES)


# Remove mpi and system libs
set(Fortran_IMPLICIT_LIBRARY_DIRS "")
set(Fortran_IMPLICIT_LIBRARY_NAMES "")
set(Fortran_IMPLICIT_LIBRARIES "")
set(Fortran_IGNORED_LIBRARIES "")
foreach (scilib ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# Whether finished with this library
  set(libprocessed FALSE)

# Ignore mpi libraries
  if (${scilib} MATCHES "^mpich")
    if (DEBUG_CMAKE)
      message("${scilib} is an MPICH library.  Ignoring.")
    endif ()
    set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
    set(libprocessed TRUE)
  elseif (${scilib} MATCHES "^open-" OR ${scilib} MATCHES "^libmpi")
    if (DEBUG_CMAKE)
      message("${scilib} is an OpenMPI library.  Ignoring.")
    endif ()
    set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
    set(libprocessed TRUE)
  endif ()

# Ignore system libraries
  if (NOT libprocessed)
    foreach (lib pthread dl nsl util rt m c)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is a system library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# Ignore communication and queue libraries
  if (NOT libprocessed)
    foreach (lib rdmacm ibverbs torque)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is an infiniband or torque library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# Ignore Hopper Cray pgi wrapper libraries
  if (NOT libprocessed)
    foreach (lib fftw3 fftw3f AtpSigHandler scicpp_pgi sci_pgi_mp mpl sma xpmem dmapp ugni pmi alpslli alpsutil udreg zceh stdmpz Cmpz pgmp nspgc pgc)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is a Hopper Cray pgi wrapper library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# Ignore additional Franklin pgi Cray wrapper libraries
  if (NOT libprocessed)
    foreach (lib scicpp stdc++ sci_quadcore_mp portals)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is a additional Franklin pgi Cray wrapper library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# Ignore additional Franklin pathscale Cray wrapper libraries
  if (NOT libprocessed)
    foreach (lib openmp eh mv mpath pscrt)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is a additional Franklin pathscale Cray wrapper library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# Ignore bgp libraries
  if (NOT libprocessed)
    foreach (lib opa dcmf.cnk dcmfcoll.cnk SPI.cna)
      if (${scilib} STREQUAL ${lib})
        if (DEBUG_CMAKE)
          message("${scilib} is a BGP compute node library.  Ignoring.")
        endif ()
        set(Fortran_IGNORED_LIBRARIES ${Fortran_IGNORED_LIBRARIES} ${scilib})
        set(libprocessed TRUE)
        break ()
      endif ()
    endforeach ()
  endif ()

# The remaining libs are added
  if (NOT libprocessed)
    set(scilibpathvar ${scilib}_LIBRARY)        # Cache variable
    # message("Looking for ${scilib}.")
    find_library(${scilibpathvar} ${scilib}
        PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}
        NO_DEFAULT_PATH)
    set(scilibpath ${${scilibpathvar}})
    if (scilibpath)
      set(Fortran_IMPLICIT_LIBRARIES ${Fortran_IMPLICIT_LIBRARIES} ${scilibpath})
      # get_filename_component(scilibname ${scilibpath} NAME_WE)
      set(Fortran_IMPLICIT_LIBRARY_NAMES ${Fortran_IMPLICIT_LIBRARY_NAMES} ${scilib})
      get_filename_component(scilibdir ${scilibpath}/.. REALPATH)
      set(Fortran_IMPLICIT_LIBRARY_DIRS ${Fortran_IMPLICIT_LIBRARY_DIRS} ${scilibdir})
    endif ()
  endif ()

endforeach ()
list(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARIES)
list(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARY_NAMES)
list(REMOVE_DUPLICATES Fortran_IMPLICIT_LIBRARY_DIRS)

# JRC, 20111203: Why are we doing this?  We can use the other variables.
if (0)
unset(Fortran_IMPLICIT_LIBFLAGS)
foreach (scilibdir ${Fortran_IMPLICIT_LIBRARY_DIRS})
  set(Fortran_IMPLICIT_LIBFLAGS
    "${Fortran_IMPLICIT_LIBFLAGS} -L${scilibdir} -Wl,-rpath,${scilibdir}")
endforeach ()
foreach (scilib ${Fortran_IMPLICIT_LIBRARY_NAMES})
  set(Fortran_IMPLICIT_LIBFLAGS
    "${Fortran_IMPLICIT_LIBFLAGS} -l${scilib}")
endforeach ()
endif ()

if (Fortran_IMPLICIT_LIBFLAGS)
  string(STRIP ${Fortran_IMPLICIT_LIBFLAGS} Fortran_IMPLICIT_LIBFLAGS)
endif ()
SciPrintString("RESULTS FOR fortran implicit libraries after removing duplicates and wrapper libraries:")
SciPrintVar(Fortran_IMPLICIT_LIBRARIES)
SciPrintVar(Fortran_IMPLICIT_LIBRARY_NAMES)
SciPrintVar(Fortran_IMPLICIT_LIBRARY_DIRS)
SciPrintVar(Fortran_IMPLICIT_LIBFLAGS)
if (DEBUG_CMAKE)
  SciPrintVar(Fortran_IGNORED_LIBRARIES)
endif ()

# Set release flags.  Assume same for now.  If different, we will
# put in the if, elseif coding.
if ("${CMAKE_BUILD_TYPE}" STREQUAL RELEASE AND OPTIMIZATION)
  set(CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
  SciPrintVar(CMAKE_Fortran_FLAGS_RELEASE)
endif ()

# Default macros for the fortran-C interface
set(FC_FUNC  "FC_FUNC(name,NAME) name ## _")
set(FC_FUNC_ "FC_FUNC_(name,NAME) name ## _")

#
# Detect name mangling convention used between Fortran and C
#
include(FortranCInterface)
FortranCInterface_HEADER(
  ${CMAKE_BINARY_DIR}/FCMangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOL_NAMESPACE "FC_"
  SYMBOLS mysub mymod:my_sub
)

file(STRINGS ${CMAKE_BINARY_DIR}/FCMangle.h CONTENTS REGEX "FC_GLOBAL\\(.*,.*\\) +(.*)")
string(REGEX MATCH "FC_GLOBAL\\(.*,.*\\) +(.*)" RESULT "${CONTENTS}")
set(FC_FUNC "FC_FUNC(name,NAME) ${CMAKE_MATCH_1}")

file(STRINGS ${CMAKE_BINARY_DIR}/FCMangle.h CONTENTS REGEX "FC_GLOBAL_\\(.*,.*\\) +(.*)")
string(REGEX MATCH "FC_GLOBAL_\\(.*,.*\\) +(.*)" RESULT "${CONTENTS}")
set(FC_FUNC_ "FC_FUNC_(name,NAME) ${CMAKE_MATCH_1}")

#
# Detect module file name
#
message(STATUS "Compiling trycompile/modulesrcfile.f90.")
execute_process(
  COMMAND ${CMAKE_Fortran_COMPILER} -c ${SCICMAKE_DIR}/trycompile/modulesrcfile.f90 -o modulesrcfile.o
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
set(SCI_FC_MODULENAME_CAPITALIZED FALSE)
set(SCI_FC_MODULE_SUFFIX)
if (EXISTS ${CMAKE_BINARY_DIR}/modulename.mod)
  set(SCI_FC_MODULE_SUFFIX mod)
  file(REMOVE ${CMAKE_BINARY_DIR}/modulename.mod)
elseif (EXISTS ${CMAKE_BINARY_DIR}/modulename.MOD)
  set(SCI_FC_MODULE_SUFFIX MOD)
  file(REMOVE ${CMAKE_BINARY_DIR}/modulename.MOD)
elseif (EXISTS ${CMAKE_BINARY_DIR}/MODULENAME.MOD)
  set(SCI_FC_MODULE_SUFFIX MOD)
  set(SCI_FC_MODULENAME_CAPITALIZED TRUE)
  file(REMOVE ${CMAKE_BINARY_DIR}/MODULENAME.MOD)
elseif (EXISTS ${CMAKE_BINARY_DIR}/MODULENAME.mod)
  set(SCI_FC_MODULE_SUFFIX MOD)
  file(REMOVE ${CMAKE_BINARY_DIR}/MODULENAME.mod)
endif ()
file(REMOVE ${CMAKE_BINARY_DIR}/modulesrcfile.o)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
SciPrintVar(Fortran_COMPILER_NAME)
SciPrintVar(SCI_FC_MODULENAME_CAPITALIZED)
SciPrintVar(SCI_FC_MODULE_SUFFIX)
set(CMAKE_Fortran_FLAGS_FULL ${CMAKE_Fortran_FLAGS_RELEASE})
SciPrintVar(CMAKE_Fortran_FLAGS_FULL)
SciPrintVar(CMAKE_Fortran_FLAGS_RELEASE)
SciPrintVar(CMAKE_Fortran_FLAGS_RELWITHDEBINFO)
SciPrintVar(CMAKE_Fortran_FLAGS_DEBUG)
SciPrintVar(CMAKE_Fortran_FLAGS)

