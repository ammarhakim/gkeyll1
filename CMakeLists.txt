# Project information
project(lucee)

# Required version
cmake_minimum_required(VERSION 2.4)

if (WIN32)
# On Windows we need to shut up security warnings
   add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif (WIN32)

# Add in path to additional modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

# Check for include files
add_definitions(-DHAVE_CONFIG_H)
INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(iostream HAVE_STD_STREAMS)
CHECK_INCLUDE_FILE_CXX(sstream HAVE_SSTREAM)

# Create configure file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Build code in loki
add_subdirectory(etc/loki/src)

# Build code in lua
add_subdirectory(etc/lua/src)

# Build code in lib
add_subdirectory(lib)

# Build all unit tests
add_subdirectory(unit)

# Add individual tests
ENABLE_TESTING()
# lcfield
add_test(lcfield unit/lcfield)

