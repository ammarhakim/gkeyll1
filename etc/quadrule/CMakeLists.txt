# Includes for this level
include_directories(${PROJECT_SOURCE_DIR}/etc/quadrule)

# Create the library
add_library(lcquadrule quadrule.c)
