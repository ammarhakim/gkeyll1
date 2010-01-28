# $Id: CMakeLists.txt 175 2008-07-13 13:57:33Z cary $

# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)


# Create the test executables
add_library(lclib 
  LcExcept.cc 
  LcKeyVal.cc
  LcKeyValTree.cc
  LcLinAlgebra.cc
  LcMatrix.cc
  LcSettableObject.cc
  LcVector.cc)
