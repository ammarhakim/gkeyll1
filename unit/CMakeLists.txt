# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/unit)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)

if (FFTW3_FOUND)
   include_directories(${Fftw3_INCLUDE_DIRS})
endif (FFTW3_FOUND)

link_directories(${PROJECT_BINARY_DIR}/slvrs)
link_directories(${PROJECT_BINARY_DIR}/datastructs)
link_directories(${PROJECT_BINARY_DIR}/grids)
link_directories(${PROJECT_BINARY_DIR}/lib)
link_directories(${PROJECT_BINARY_DIR}/etc/loki/src)

# add the Accelerate framework on Mac OS X
IF(APPLE)
  FIND_LIBRARY(ACCELERATE_LIB Accelerate)
  MARK_AS_ADVANCED(ACCELERATE_LIB)
  SET(EXTRA_LIBS ${ACCELERATE_LIB})
ELSE (APPLE)
  IF (LAPACK_FOUND)
    SET(EXTRA_LIBS ${Lapack_LIBRARIES})
  ENDIF ()
ENDIF (APPLE)

# Create test executables

if (ENABLE_PARALLEL)
  set (Hdf5_Mpi_LIBRARIES ${TxBase_LIBRARIES} ${Hdf5_LIBRARIES} ${MPI_STLIBS} ${Z_LIBRARIES})
else ()
  set (Hdf5_Mpi_LIBRARIES ${TxBase_LIBRARIES} ${Hdf5_LIBRARIES})
endif ()

# lcarray
add_executable(lcarray lcarray.cxx)
target_link_libraries(lcarray lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES})

# lccolmajorindexer
add_executable(lccolmajorindexer lccolmajorindexer.cxx)
target_link_libraries(lccolmajorindexer lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcrowmajorindexer
add_executable(lcrowmajorindexer lcrowmajorindexer.cxx)
target_link_libraries(lcrowmajorindexer lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcdirsequencer
add_executable(lcdirsequencer lcdirsequencer.cxx)
target_link_libraries(lcdirsequencer lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcfixedvector
add_executable(lcfixedvector lcfixedvector.cxx)
target_link_libraries(lcfixedvector lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcvector
add_executable(lcvector lcvector.cxx)
target_link_libraries(lcvector lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES})

# lcmatrix
add_executable(lcmatrix lcmatrix.cxx)
target_link_libraries(lcmatrix lclib lcloki ${Lua_STLIBS} ${EXTRA_LIBS} ${Hdf5_Mpi_LIBRARIES})

# lccmdlineargs
add_executable(lccmdlineargs lccmdlineargs.cxx)
target_link_libraries(lccmdlineargs lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcregion
add_executable(lcregion lcregion.cxx)
target_link_libraries(lcregion lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcfield
add_executable(lcfield lcfield.cxx)
target_link_libraries(lcfield
  lcdatastructs lcgrids lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcrowmajorsequencer
add_executable(lcrowmajorsequencer lcrowmajorsequencer.cxx)
target_link_libraries(lcrowmajorsequencer lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lccolmajorsequencer
add_executable(lccolmajorsequencer lccolmajorsequencer.cxx)
target_link_libraries(lccolmajorsequencer lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcluatable
add_executable(lcluatable lcluatable.cxx)
target_link_libraries(lcluatable lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcluamodule
add_executable(lcluamodule lcluamodule.cxx)
target_link_libraries(lcluamodule lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lclogger
add_executable(lclogger lclogger.cxx)
target_link_libraries(lclogger lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lchdf5io
add_executable(lchdf5io lchdf5io.cxx)
target_link_libraries(lchdf5io lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES})

# lcvec3
add_executable(lcvec3 lcvec3.cxx)
target_link_libraries(lcvec3 lclib ${Hdf5_Mpi_LIBRARIES})

# lcrectcartgrid
add_executable(lcrectcartgrid lcrectcartgrid.cxx)
target_link_libraries(lcrectcartgrid lcgrids lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lclincombiner
add_executable(lclincombiner lclincombiner.cxx)
target_link_libraries(lclincombiner lcslvrs lcdatastructs 
  lcgrids lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES} ${Gsl_LIBRARIES})

# lcluaobj
add_executable(lcluacobj lcluacobj.cxx)
target_link_libraries(lcluacobj lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcfieldptr
add_executable(lcfieldptr lcfieldptr.cxx)
target_link_libraries(lcfieldptr lcdatastructs lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcgslelliptic
add_executable(lcgslelliptic lcgslelliptic.cxx)
target_link_libraries(lcgslelliptic ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcrectcoordsys
add_executable(lcrectcoordsys lcrectcoordsys.cxx)
target_link_libraries(lcrectcoordsys lcdatastructs lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcvalarraycmp
add_executable(lcvalarraycmp lcvalarraycmp.cxx)
target_link_libraries(lcvalarraycmp lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES})

# lcgeoprim
add_executable(lcgeomprim lcgeomprim.cxx)
target_link_libraries(lcgeomprim lcgrids lclib lcloki ${Lua_STLIBS} ${Hdf5_Mpi_LIBRARIES})

# lcstructgridfield
add_executable(lcstructgridfield lcstructgridfield.cxx)
target_link_libraries(lcstructgridfield lcdatastructs lcgrids lclib lcloki ${Lua_STLIBS} 
  ${Hdf5_Mpi_LIBRARIES} ${Gsl_LIBRARIES})

# lcgridode
add_executable(lcgridode lcgridode.cxx)
target_link_libraries(lcgridode lcslvrs lcdatastructs lcgrids lclib lcloki ${Lua_STLIBS}
  ${EXTRA_LIBS} ${Hdf5_Mpi_LIBRARIES} ${Gsl_LIBRARIES})

if (FFTW3_FOUND)
  add_executable(lcfftwtest lcfftwtest.cxx)
  target_link_libraries(lcfftwtest ${Fftw3_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})
endif (FFTW3_FOUND)

# lcdecompregion
add_executable(lcdecompregion lcdecompregion.cxx)
target_link_libraries(lcdecompregion lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcmultiregion
add_executable(lcmultiregion lcmultiregion.cxx)
target_link_libraries(lcmultiregion lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcquadraturerule
add_executable(lcquadraturerule lcquadraturerule.cxx)
target_link_libraries(lcquadraturerule lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})

# lcdynvector
add_executable(lcdynvector lcdynvector.cxx)
target_link_libraries(lcdynvector 
  lcdatastructs lclib lcloki ${Lua_STLIBS} ${Gsl_LIBRARIES} ${Hdf5_Mpi_LIBRARIES})
