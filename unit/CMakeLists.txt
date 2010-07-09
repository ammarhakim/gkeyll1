# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/unit)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

link_directories(${PROJECT_BINARY_DIR}/slvrs)
link_directories(${PROJECT_BINARY_DIR}/grids)
link_directories(${PROJECT_BINARY_DIR}/lib)
link_directories(${PROJECT_BINARY_DIR}/etc/loki/src)
link_directories(${PROJECT_BINARY_DIR}/etc/lua/src)

# add the Accelerate framework on Mac OS X
IF(APPLE)
  FIND_LIBRARY(ACCELERATE_LIB Accelerate)
  MARK_AS_ADVANCED(ACCELERATE_LIB)
  SET(EXTRA_LIBS ${ACCELERATE_LIB})
ENDIF (APPLE)

# Create test executables

# lcarray
add_executable(lcarray lcarray.cxx)
target_link_libraries(lcarray lclib lcloki lclua)

# lccolmajorindexer
add_executable(lccolmajorindexer lccolmajorindexer.cxx)
target_link_libraries(lccolmajorindexer lclib lcloki lclua)

# lcrowmajorindexer
add_executable(lcrowmajorindexer lcrowmajorindexer.cxx)
target_link_libraries(lcrowmajorindexer lclib lcloki lclua)

# lcfixedvector
add_executable(lcfixedvector lcfixedvector.cxx)
target_link_libraries(lcfixedvector lclib lcloki lclua)

# lcvector
add_executable(lcvector lcvector.cxx)
target_link_libraries(lcvector lclib lcloki lclua)

# lcmatrix
add_executable(lcmatrix lcmatrix.cxx)
target_link_libraries(lcmatrix lclib lcloki lclua ${EXTRA_LIBS})

# lccmdlineargs
add_executable(lccmdlineargs lccmdlineargs.cxx)
target_link_libraries(lccmdlineargs lclib lcloki lclua)

# lcregion
add_executable(lcregion lcregion.cxx)
target_link_libraries(lcregion lclib lcloki lclua)

# lcfield
add_executable(lcfield lcfield.cxx)
target_link_libraries(lcfield lclib lcloki lclua)

# lcrowmajorsequencer
add_executable(lcrowmajorsequencer lcrowmajorsequencer.cxx)
target_link_libraries(lcrowmajorsequencer lclib lcloki lclua)

# lccolmajorsequencer
add_executable(lccolmajorsequencer lccolmajorsequencer.cxx)
target_link_libraries(lccolmajorsequencer lclib lcloki lclua)

# lcobjcreator
add_executable(lcobjcreator lcobjcreator.cxx)
target_link_libraries(lcobjcreator lclib lcloki lclua)

# lcluatable
add_executable(lcluatable lcluatable.cxx)
target_link_libraries(lcluatable lclib lcloki lclua)

# lcluamodule
add_executable(lcluamodule lcluamodule.cxx)
target_link_libraries(lcluamodule lclib lcloki lclua)

# lclogger
add_executable(lclogger lclogger.cxx)
target_link_libraries(lclogger lclib lcloki lclua)

# lchdf5io
add_executable(lchdf5io lchdf5io.cxx)
target_link_libraries(lchdf5io lclib lcloki lclua ${HDF5_LIBRARIES})

# lcvec3
add_executable(lcvec3 lcvec3.cxx)
target_link_libraries(lcvec3 lclib)

# lcrectcartgrid
add_executable(lcrectcartgrid lcrectcartgrid.cxx)
target_link_libraries(lcrectcartgrid lcgrids lclib)

# lcluafunction
add_executable(lcluafunction lcluafunction.cxx)
target_link_libraries(lcluafunction lclib lcloki lclua)

# lclincombiner
add_executable(lclincombiner lclincombiner.cxx)
target_link_libraries(lclincombiner lcslvrs lcgrids lclib lcloki lclua)
