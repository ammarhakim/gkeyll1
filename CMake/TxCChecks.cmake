######################################################################
#
# TxCChecks: check various C capabilities
#
# $Id: TxCChecks.cmake 182 2010-12-29 17:13:03Z cary $
#
# Copyright 2010 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

# Determine the compiler id
SET(C_COMPILER ${CMAKE_C_COMPILER})
SET(C_COMPILER_ID ${CMAKE_C_COMPILER_ID})
MESSAGE(STATUS "C_COMPILER_ID = ${C_COMPILER_ID}.")

# Check whether time and sys/time can both be included
INCLUDE(CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES(
"#include <sys/time.h>\n#include <time.h>\nint main(int argc, char** argv) {return 0;}"
TIME_WITH_SYS_TIME
)
IF(TIME_WITH_SYS_TIME)
  IF(DEBUG_CMAKE)
    MESSAGE("time.h and sys/time.h are compatible.")
  ENDIF()
  SET(TIME_WITH_SYS_TIME 1 CACHE BOOL "Whether time and sys/time are compatible")
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("time.h and sys/time.h are NOT compatible.")
  ENDIF()
ENDIF()

# Check whether struct tm is in sys/time
INCLUDE(CheckCSourceCompiles)
CHECK_C_SOURCE_COMPILES(
"#include <sys/time.h>\nint main(int argc, char** argv) {struct tm tm;int *p = &tm.tm_sec;return !p;}"
TM_IN_SYS_TIME
)
IF(TM_IN_SYS_TIME)
  IF(DEBUG_CMAKE)
    MESSAGE("struct tm is in time.h.")
  ENDIF()
  SET(TM_IN_SYS_TIME 1 CACHE BOOL "Whether struct tm is in sys/time.")
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("struct tm is NOT in time.h.")
  ENDIF()
ENDIF()

# Check variable sizes
CHECK_C_SOURCE_COMPILES(
"#include <stdio.h>\nvoid f(unsigned int i){}\nvoid f(size_t i){}\nint main(int argc, char** argv) {return 0;}"
UINT_IS_NOT_SIZE_T
)
IF(UINT_IS_NOT_SIZE_T)
  IF(DEBUG_CMAKE)
    MESSAGE("uint and size_t are not the same size.")
  ENDIF()
  SET(UINT_IS_SIZE_T 1 CACHE BOOL "Whether uint is the same as size_t")
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("uint and size_t are the same size.")
  ENDIF()
ENDIF()

CHECK_C_SOURCE_COMPILES(
"#include <stdio.h>\nvoid f(int i){}\nvoid f(ssize_t i){}\nint main(int argc, char** argv) {return 0;}"
INT_IS_NOT_SSIZE_T
)
IF(INT_IS_NOT_SSIZE_T)
  IF(DEBUG_CMAKE)
    MESSAGE("int and ssize_t are not the same size.")
  ENDIF()
  SET(INT_IS_SSIZE_T 1 CACHE BOOL "Whether int is the same as ssize_t")
ELSE()
  IF(DEBUG_CMAKE)
    MESSAGE("int and ssize_t are the same size.")
  ENDIF()
ENDIF()

# Get math into C for Windows
IF(WIN32)
  SET(_USE_MATH_DEFINES 1 CACHE BOOL "To bring in math defines on Windows.")
ENDIF()

# Determine RELEASE build flags based on optimization variable.
# Correspondence
#          autotools                      cmake
# --with-optimization=default   -DCMAKE_BUILD_TYPE:STRING=RELWITHDEBINFO
# --with-optimization=debug     -DCMAKE_BUILD_TYPE:STRING=DEBUG
# Skip for CMAKE_BUILD_TYPE:STRING = DEBUG or RELWITHDEBINFO.
# Debug build should be obtained by setting -DCMAKE_BUILD_TYPE:STRING=DEBUG

# For default, nothing is set
IF(${CMAKE_BUILD_TYPE} STREQUAL RELEASE AND OPTIMIZATION)

# Determine the processor and sse capabilities
  IF(EXISTS /proc/cpuinfo)
    MESSAGE(STATUS "Working on LINUX")
    FILE(READ /proc/cpuinfo cpuinfo)
    EXECUTE_PROCESS(COMMAND cat /proc/cpuinfo
      COMMAND grep "model name"
      COMMAND head -1
      OUTPUT_VARIABLE TXC_CPU
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    STRING(REGEX REPLACE "^.*: " "" TXC_CPU ${TXC_CPU})
    MESSAGE(STATUS "CPU is ${TXC_CPU}")
    EXECUTE_PROCESS(COMMAND cat /proc/cpuinfo
      COMMAND grep "flags"
      COMMAND head -1
      OUTPUT_VARIABLE CPU_CAPABILITIES
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ELSEIF(APPLE)
    EXECUTE_PROCESS(COMMAND sysctl -a machdep.cpu.brand_string
      OUTPUT_VARIABLE TXC_CPU
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    STRING(REGEX REPLACE "^.*: " "" TXC_CPU ${TXC_CPU})
    MESSAGE(STATUS "CPU is ${TXC_CPU}")
    EXECUTE_PROCESS(COMMAND sysctl -a machdep.cpu.features
      COMMAND head -1
      OUTPUT_VARIABLE CPU_CAPABILITIES
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    STRING(REGEX REPLACE "SSE" "sse" CPU_CAPABILITIES ${CPU_CAPABILITIES})
  ENDIF()
  IF(CPU_CAPABILITIES)
    STRING(REGEX REPLACE "^.*: *" "" CPU_CAPABILITIES ${CPU_CAPABILITIES})
    SEPARATE_ARGUMENTS(CPU_CAPABILITIES)
    MESSAGE(STATUS "CPU capabilities are ${CPU_CAPABILITIES}")
    FOREACH(cap ${CPU_CAPABILITIES})
      # MESSAGE("Examining ${cap}")
      IF(${cap} MATCHES "^sse.*$")
        LIST(APPEND SSE_CAPABILITIES ${cap})
      ENDIF()
    ENDFOREACH()
    LIST(SORT SSE_CAPABILITIES)
    LIST(REVERSE SSE_CAPABILITIES)
    MESSAGE(STATUS "SSE capabilities are ${SSE_CAPABILITIES}")
    LIST(GET SSE_CAPABILITIES 0 SSE_CAPABILITY)
    MESSAGE(STATUS "SSE capability is ${SSE_CAPABILITY}")
  ENDIF()

  IF(${C_COMPILER_ID} STREQUAL GNU)
    IF(${OPTIMIZATION} STREQUAL full)
      SET(CMAKE_C_FLAGS_RELEASE "-O3 -ffast-math")
      IF(SSE_CAPABILITY)
        SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -m${SSE_CAPABILITY}")
      ENDIF()
# Apple automatically sets mtune
      # IF(LINUX)
        IF("${TXC_CPU}" MATCHES ".*Athlon.*MP.*")
          SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=athlon-mp")
        ELSEIF(${TXC_CPU} MATCHES ".*Athlon.*XP.*")
          SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=athlon-xp")
        ELSEIF(${TXC_CPU} MATCHES ".*Athlon.*")
          SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=athlon")
        ELSEIF(${TXC_CPU} MATCHES ".*Opteron.*")
          SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mtune=amdfam10")
        ENDIF()
      # ENDIF()
    ENDIF()

  ELSEIF(${C_COMPILER_ID} STREQUAL Intel)

  ELSEIF(${C_COMPILER_ID} STREQUAL MSVC)

  ELSEIF(${C_COMPILER_ID} STREQUAL PathScale)

  ELSEIF(${C_COMPILER_ID} STREQUAL PGI)

  ELSEIF(${C_COMPILER_ID} STREQUAL XL)

  ELSE()
    MESSAGE(STATUS "Flags not known for ${C_COMPILER_ID}")
  ENDIF()

  MESSAGE(STATUS "CMAKE_C_FLAGS_RELEASE = '${CMAKE_C_FLAGS_RELEASE}'")

ENDIF()

