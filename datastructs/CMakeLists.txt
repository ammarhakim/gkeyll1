# Includes for this level
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/lib)
include_directories(${PROJECT_SOURCE_DIR}/grids)
include_directories(${PROJECT_SOURCE_DIR}/datastructs)
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)
include_directories(${PROJECT_SOURCE_DIR}/etc/lua/src)


# Create the library
add_library(lcdatastructs
  LcConstFieldPtr.cpp
  LcDataStructRegistry.cpp
  LcDynVector.cpp
  LcField.cpp
  LcFieldPtr.cpp
  LcStructGridField.cpp
)

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY ./
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
  )
  
  install(TARGETS lcdatastructs
    ARCHIVE DESTINATION lib
  )
endif ()
