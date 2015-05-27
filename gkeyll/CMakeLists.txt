# Build gkeyll executable

# Project include directories
INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/gkeyll
  ${CMAKE_SOURCE_DIR}/rte
  ${CMAKE_SOURCE_DIR}/proto
  ${CMAKE_SOURCE_DIR}/pic
  ${CMAKE_SOURCE_DIR}/slvrs
  ${CMAKE_SOURCE_DIR}/lib
  ${CMAKE_SOURCE_DIR}/grids
  ${CMAKE_SOURCE_DIR}/datastructs
  ${CMAKE_SOURCE_DIR}/etc/quadrule
  ${CMAKE_SOURCE_DIR}/etc/loki/include
  ${CMAKE_BINARY_DIR}
)

set(lucee_LIBRARIES
  lcrte
  lcproto
  lcslvrs
  lcpic
  lcdatastructs
  lcgrids
  lclib
  lcquadrule
  lcloki
)

#
# External includes and libraries
#
unset(FC_EXTERNAL_INCDIRS)
unset(FC_EXTERNAL_LIBFLAGS)
unset(FC_EXTERNAL_LIBDIRS)
unset(FacetsSt_EXT_LIBS)

if (FFTW3_FOUND)
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${Fftw3_INCLUDE_DIRS})
  set(FC_EXTERNAL_FLAGS "${FC_EXTERNAL_FLAGS} ${Fftw3_FLAGS}")
  set(FC_EXTERNAL_LIBS ${FC_EXTERNAL_LIBS} ${Fftw3_LIBRARIES})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${Fftw3_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Fftw3_STLIBS})
endif ()
 
if (GSL_FOUND)
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${Gsl_INCLUDE_DIRS})
  set(FC_EXTERNAL_FLAGS "${FC_EXTERNAL_FLAGS} ${Gsl_FLAGS}")
  set(FC_EXTERNAL_LIBS ${FC_EXTERNAL_LIBS} ${Gsl_LIBRARIES})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${Gsl_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Gsl_STLIBS})
endif ()

if (PETSC_FOUND)
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${Petsc_INCLUDE_DIRS})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${Petsc_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Petsc_STLIBS})
  set(FacetsDl_EXT_LIBS ${FacetsDl_EXT_LIBS} ${Petsc_STLIBS})
endif ()

# clapack only on non-Mac platforms
if (LAPACK_FOUND AND (NOT APPLE))
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${Lapack_INCLUDE_DIRS})
  set(FC_EXTERNAL_FLAGS "${FC_EXTERNAL_FLAGS} ${Lapack_FLAGS}")
  set(FC_EXTERNAL_LIBS ${FC_EXTERNAL_LIBS} ${Lapack_LIBRARIES})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${Lapack_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Lapack_STLIBS})
endif ()

if (MOAB_FOUND)
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${Moab_INCLUDE_DIRS})
  set(FC_EXTERNAL_FLAGS "${FC_EXTERNAL_FLAGS} ${Moab_FLAGS}")
  set(FC_EXTERNAL_LIBS ${FC_EXTERNAL_LIBS} ${Moab_LIBRARIES})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${Moab_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Moab_STLIBS})
endif ()

# Required base directories
set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS}
  ${Lua_INCLUDE_DIRS}
)
set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS}
  ${Lua_LIBRARY_DIRS}
)
set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS}
  ${Lua_STLIBS}
)

set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS}
  ${TxBase_INCLUDE_DIRS} ${Hdf5_INCLUDE_DIRS}
)
set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS}
  ${TxBase_LIBRARY_DIRS} ${Hdf5_LIBRARY_DIRS}
)
set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS}
  ${TxBase_STLIBS} ${Hdf5_STLIBS}
)

if (DEBUG_CMAKE)
  message("nautilus/CMakeLists.txt: Adding ${TxBase_STLIBS} ${Hdf5_STLIBS} to FacetsSt_EXT_LIBS.")
endif ()

# explicitly add in blas framework if on apple
if (APPLE)
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} "-framework Accelerate")
endif ()

# Add in fortran libs
set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${Fortran_IMPLICIT_LIBRARIES})

if (ENABLE_PARALLEL AND MPI_INCLUDE_PATH AND NOT TX_HAVE_MPICXX_COMPILER_WRAPPER)
  set(FC_EXTERNAL_INCDIRS ${FC_EXTERNAL_INCDIRS} ${MPI_INCLUDE_DIRS})
  set(FC_EXTERNAL_LIBS ${FC_EXTERNAL_LIBS} ${MPI_LIBRARIES})
  set(FC_EXTERNAL_LIBDIRS ${FC_EXTERNAL_LIBDIRS} ${MPI_LIBRARY_DIRS})
  set(FacetsSt_EXT_LIBS ${FacetsSt_EXT_LIBS} ${MPI_STLIBS})
endif ()

SciPrintString("")
SciPrintString("-------- FACETS external and system libraries --------")
string(STRIP "${FC_EXTERNAL_LIBFLAGS}" FC_EXTERNAL_LIBFLAGS)
SciPrintVar(FC_EXTERNAL_LIBFLAGS)
SciPrintVar(FC_EXTERNAL_LIBDIRS)
SciPrintVar(FacetsSt_EXT_LIBS)

# Add in all include directories
include_directories(${FC_EXTERNAL_INCDIRS})

# add all external library dirs to link
link_directories(${FC_EXTERNAL_LIBDIRS})

if (ENABLE_PARALLEL AND MPI_INCLUDE_PATH AND NOT TX_HAVE_MPICXX_COMPILER_WRAPPER)
  include_directories("${MPI_INCLUDE_PATH}")
endif ()

# System libraries that we can assume always present
# Python needs these libraries on Hopper, but they have to be at the end.
set(FC_SYS_LIBS)
if (WIN32)
  set(FC_SYS_LIBS ws2_32 winmm)
else ()
  # if (PYTHONLIBS_FOUND)
# These are not needed if one is using the openmpi compiler wrappers.
# But on a Cray or elsewhere these must be added.
    set(FC_SYS_LIBS ${FC_SYS_LIBS} pthread util dl)
  # endif ()
# pgi-11.2 is getting some implementations from the system runtime library
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL PGI OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL XL)
    set(FC_SYS_LIBS ${FC_SYS_LIBS} rt)
  endif ()
endif ()
SciPrintVar(FC_SYS_LIBS)
SciGetStaticLibs("${FC_SYS_LIBS}" FC_SYS_STLIBS)
SciPrintVar(FC_SYS_STLIBS)

#
# Executables
#
if (ENABLE_PARALLEL)
  set(LUCEE_STATIC gkeyll)
else ()
  set(LUCEE_STATIC gkeyllser)
endif ()

# Static link
add_executable(${LUCEE_STATIC} gkeyll.cxx LcRegisterModules.cpp)
set_target_properties(${LUCEE_STATIC} PROPERTIES
  LINK_FLAGS "${FC_EXTERNAL_FLAGS}"
)
target_link_libraries(${LUCEE_STATIC}  ${lucee_LIBRARIES})
target_link_libraries(${LUCEE_STATIC} ${FacetsSt_EXT_LIBS} )
#${Moab_LIBRARIES}

message("The xlibraries are : " ${Moab_LIBRARIES})

if (USE_STATIC_SYSLIBS)
  message(STATUS "Using static system libraries: ${FC_SYS_STLIBS}")
  target_link_libraries(${LUCEE_STATIC} ${FC_SYS_STLIBS})
# This gets the link line to end with -Wl,-Bstatic, so that subsequent
# libraries linked by the Cray compiler wrappers are linked static.
  target_link_libraries(${LUCEE_STATIC} -Wl,-Bstatic)
else ()
  target_link_libraries(${LUCEE_STATIC} ${FC_SYS_LIBS})
  message(STATUS "Not using static system libraries: ${FC_SYS_LIBS}")
  if (WIN32)
  elseif (APPLE)
  else ()
    set_target_properties(${LUCEE_STATIC} PROPERTIES INSTALL_RPATH "\$ORIGIN:\$ORIGIN/../lib")
  endif ()
endif ()

install(TARGETS ${LUCEE_STATIC}
  DESTINATION bin
  PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
              GROUP_READ GROUP_EXECUTE ${TX_GROUP_WRITE}
)

