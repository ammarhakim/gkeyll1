# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)


# Create the library
add_library(lcgrids
  LcCellElem.cpp
  LcEdgeElem.cpp
  LcFaceElem.cpp
  LcGmvGridCreator.cpp
  LcGridFactory.cpp
  LcGridGeometry.cpp
  LcGridIfc.cpp
  LcGridRegistry.cpp
  LcMappedCartGrid.cpp
  LcMappedCartGridFactory.cpp
  LcRectCartGrid.cpp
  LcRectCartGridFactory.cpp
  LcStructuredGridBase.cpp
  LcUnstructConnectivity.cpp
  LcUnstructGrid.cpp
  LcUnstructGridCreator.cpp
  LcVertexElem.cpp
)
