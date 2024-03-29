# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/libb64/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

# Create the library
add_library(lclib 
  LcAlignedRectCoordSys.cpp
  LcBasicObj.cpp
  LcCartGeneralDecompRegionCalc.cpp
  LcCartProdDecompRegionCalc.cpp
  LcCmdLineArgs.cpp
  LcDataStructIfc.cpp
  LcDecompRegion.cpp
  LcDecompRegionCalcIfc.cpp
  LcDecompRegionCalcRegistry.cpp
  LcExcept.cpp 
  LcFileHandler.cpp
  LcGaussianQuadRule.cpp
  LcLibRegistry.cpp
  LcLinAlgebra.cpp
  LcLogRecordHandler.cpp
  LcLogStream.cpp
  LcLogStreamStrm.cpp
  LcLogger.cpp
  LcLuaFuncMap.cpp
  LcLuaMathLib.cpp
  LcLuaObjTypeId.cpp
  LcLuaState.cpp
  LcLuaTable.cpp
  LcLuceeMod.cpp
  LcMathLib.cpp
  LcMatrix.cpp
  LcMultiRegion.cpp
  LcQuadratureRule.cpp
  LcRectCoordSys.cpp
  LcSolverIfc.cpp
  LcStreamHandler.cpp
  LcSubCartProdDecompRegionCalc.cpp
  LcVec3.cpp
  LcVector.cpp
)

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY ./
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
  )
  
  install(TARGETS lclib
    ARCHIVE DESTINATION lib
  )
endif ()
