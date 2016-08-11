# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/loki/include)

# Create the library
add_library(lcloki
  OrderedStatic.cpp
  SafeFormat.cpp
  Singleton.cpp
  SmallObj.cpp
  SmartPtr.cpp
  StrongPtr.cpp)

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY
    ${PROJECT_SOURCE_DIR}/etc/loki/include/loki
    DESTINATION include
  )
  
  install(TARGETS lcloki
    ARCHIVE DESTINATION lib
  )
endif ()
