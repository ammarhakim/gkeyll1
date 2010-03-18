# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

# Create the library
add_library(lclib 
  LcCmdLineArgs.cpp
  LcConstFieldPtr.cpp
  LcExcept.cpp 
  LcField.cpp
  LcFieldPtr.cpp
  LcFileHandler.cpp
  LcHdf5Io.cpp
  LcHdf5IoTmpl.cpp
  LcIoBase.cpp
  LcLinAlgebra.cpp
  LcLogRecordHandler.cpp
  LcLogStream.cpp
  LcLogStreamStrm.cpp
  LcLogger.cpp
  LcLuaState.cpp
  LcLuaTable.cpp
  LcMathLib.cpp
  LcMatrix.cpp
  LcStreamHandler.cpp
  LcVector.cpp
)
