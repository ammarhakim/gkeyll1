# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)


# Create the library
add_library(lclib 
  LcCmdLineArgs.cpp
  LcExcept.cpp 
  LcField.cpp
  LcKeyVal.cpp
  LcKeyValTree.cpp
  LcLinAlgebra.cpp
  LcLuaSaveRstrStack.cpp
  LcLuaState.cpp
  LcMatrix.cpp
  LcSettableObject.cpp
  LcVector.cpp
)
