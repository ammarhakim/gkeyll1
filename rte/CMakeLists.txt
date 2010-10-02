# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lucee)
include_directories(${PROJECT_SOURCE_DIR}/rte)
include_directories(${PROJECT_SOURCE_DIR}/slvrs)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)

# Create the library
add_library(lcrte
  LcHgPhaseFunction.cpp
  LcPlCoeffsPhaseFunction.cpp
  LcRteHomogeneousSlab.cpp
  LcRtePhaseFunction.cpp
  LcRteRegistry.cpp
)


