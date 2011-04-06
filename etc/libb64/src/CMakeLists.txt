# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/libb64/include)

# Create the library
add_library(lclibb64
  cdecode.c
  cencode.c
)
