# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/quadrule)

# Create the library
add_library(lcquadrule quadrule.c)

if (NOT ${ENABLE_PARALLEL})
  install(DIRECTORY ./
    DESTINATION include
    FILES_MATCHING PATTERN "*.h"
  )
  
  install(TARGETS lcquadrule
    ARCHIVE DESTINATION lib
  )
endif ()
