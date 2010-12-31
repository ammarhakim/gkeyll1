######################################################################
#
# TxCxxChecks: check various C++ capabilities
#
# $Id: TxCxxChecks.cmake 181 2010-12-29 15:42:45Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

INCLUDE(CheckIncludeFileCXX)

CHECK_INCLUDE_FILE_CXX(iostream HAVE_STD_STREAMS)
CHECK_INCLUDE_FILE_CXX(sstream HAVE_SSTREAM)

# Determine compiler version
INCLUDE(${PROJECT_SOURCE_DIR}/CMake/TxCxxFindVersion.cmake)
IF(CXX_VERSION)
  MESSAGE(STATUS "CXX_VERSION = ${CXX_VERSION}")
ELSE(CXX_VERSION)
  MESSAGE(FATAL_ERROR "Could not determine compiler version.")
ENDIF(CXX_VERSION)

# Set the lib subdir from the Compiler ID and version
IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL GNU)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pipe")
  STRING(SUBSTRING ${CXX_VERSION} 0 1 CXX_MAJOR_VERSION)
  SET(CXX_COMP_LIB_SUBDIR gcc${CXX_MAJOR_VERSION})
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL PathScale)
  STRING(SUBSTRING ${CXX_VERSION} 0 1 CXX_MAJOR_VERSION)
  SET(CXX_COMP_LIB_SUBDIR path${CXX_MAJOR_VERSION})
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL PGI)
  # STRING(REGEX REPLACE "\\..*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  # STRING(REGEX MATCH "^[0-9]+" CXX_MAJOR_VERSION ${CXX_VERSION})
  STRING(REGEX REPLACE "\\.[0-9]+-.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  # SET(CXX_COMP_LIB_SUBDIR pgCC)  # OLD
  SET(CXX_COMP_LIB_SUBDIR pgi${CXX_MAJOR_VERSION})
ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL XL)
# This should be the basename of the compiler
  STRING(REGEX REPLACE "\\.[0-9]+.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  STRING(REGEX REPLACE "^0+" "" CXX_MAJOR_VERSION ${CXX_MAJOR_VERSION})
  GET_FILENAME_COMPONENT(REL_CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} NAME)
# Since we install ben builds in a completely different directory, can
# use same name for CXX_COMP_LIB_SUBDIR
  # IF(${REL_CMAKE_CXX_COMPILER} MATCHES "bg.*_r$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "/.*/bg.*_r$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "mp.*_r$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "/.*/mp.*_r$")
  #   SET(CXX_COMP_LIB_SUBDIR bgxlC_r${CXX_MAJOR_VERSION})
  # ELSEIF(${REL_CMAKE_CXX_COMPILER} MATCHES "bg.*$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "/.*/bg.*$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "mp.*$" OR ${REL_CMAKE_CXX_COMPILER} MATCHES "/.*/mp.*$")
  #   SET(CXX_COMP_LIB_SUBDIR bgxlC${CXX_MAJOR_VERSION})
  # ELSEIF(${REL_CMAKE_CXX_COMPILER} MATCHES ".*_r$")
  IF(${REL_CMAKE_CXX_COMPILER} MATCHES ".*_r$")
    SET(CXX_COMP_LIB_SUBDIR xlC_r${CXX_MAJOR_VERSION})
  ELSE()
    SET(CXX_COMP_LIB_SUBDIR xlC${CXX_MAJOR_VERSION})
  ENDIF()
  # SET(SEPARATE_INSTANTIATIONS 1)
  SET(SEPARATE_INSTANTIATIONS 1 CACHE BOOL "Whether to separate instantiations -- for correct compilation on xl")
ENDIF()
MESSAGE(STATUS "CXX_COMP_LIB_SUBDIR = ${CXX_COMP_LIB_SUBDIR}")

# See whether generally declared statics work
TRY_COMPILE(HAVE_GENERALLY_DECLARED_STATICS ${PROJECT_BINARY_DIR}/CMake
  ${PROJECT_SOURCE_DIR}/CMake/gendeclstatics.cpp)
IF(HAVE_GENERALLY_DECLARED_STATICS)
  IF(DEBUG_CMAKE)
    MESSAGE("${PROJECT_SOURCE_DIR}/CMake/gendeclstatics.cpp compiled.")
  ENDIF()
  SET(HAVE_GENERALLY_DECLARED_STATICS 1 CACHE BOOL "Whether the C++ compiler allows generally declared templated static variables")
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("${PROJECT_SOURCE_DIR}/CMake/gendeclstatics.cpp did not compile.")
  ENDIF()
ENDIF()

# See whether std::abs<double> known.
TRY_COMPILE(HAVE_STD_ABS_DOUBLE ${PROJECT_BINARY_DIR}/CMake
  ${PROJECT_SOURCE_DIR}/CMake/stdabsdbl.cpp)
IF(HAVE_STD_ABS_DOUBLE)
  IF(DEBUG_CMAKE)
    MESSAGE("${PROJECT_SOURCE_DIR}/CMake/stdabsdbl.cpp compiled.")
  ENDIF()
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("${PROJECT_SOURCE_DIR}/CMake/stdabsdbl.cpp did not compile.")
  ENDIF()
  SET(NOT_HAVE_STD_ABS_DOUBLE 1 CACHE BOOL "Define when the C++ compiler does not understand std::abs with double arg")
ENDIF()

# Set release flags.  Assume same for now.  If different, we will
# put in the if, elseif coding.
IF(${CMAKE_BUILD_TYPE} STREQUAL RELEASE AND OPTIMIZATION)
  SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
  MESSAGE(STATUS "CMAKE_CXX_FLAGS_RELEASE = '${CMAKE_CXX_FLAGS_RELEASE}'")
ENDIF()

