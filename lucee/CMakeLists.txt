# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lucee)
include_directories(${PROJECT_SOURCE_DIR}/rte)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

link_directories(${PROJECT_BINARY_DIR}/lib)
link_directories(${PROJECT_BINARY_DIR}/rte)
link_directories(${PROJECT_BINARY_DIR}/slvr)
link_directories(${PROJECT_BINARY_DIR}/etc/loki/src)
link_directories(${PROJECT_BINARY_DIR}/etc/lua/src)

# add the Accelerate framework on Mac OS X
IF(APPLE)
  FIND_LIBRARY(ACCELERATE_LIB Accelerate)
  MARK_AS_ADVANCED(ACCELERATE_LIB)
  SET(EXTRA_LIBS ${ACCELERATE_LIB})
ENDIF (APPLE)

# Create lucee executable
add_executable(lucee lucee.cxx LcRegisterModules.cpp)
target_link_libraries(lucee lcrte lcslvrs lclib lcloki lclua ${GSL_CONFIG_LIBS})
