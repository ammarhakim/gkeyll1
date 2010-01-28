# $Id: CMakeLists.txt 289 2009-11-17 19:53:28Z mkoch $

# Project information
project(lucee)

# Required version
cmake_minimum_required(VERSION 2.4)

# Add in the path to additional modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/Modules/")

# Check for include files
add_definitions(-DHAVE_CONFIG_H)
INCLUDE(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX(iostream HAVE_STD_STREAMS)
CHECK_INCLUDE_FILE_CXX(sstream HAVE_SSTREAM)

# Create the configure file
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config-cmake.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# Directory structure
add_subdirectory(lib)

