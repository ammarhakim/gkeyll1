######################################################################
#
# SciCxxChecks: check various C++ capabilities
#
# $Id: SciCxxChecks.cmake 656 2014-10-25 14:14:53Z jrobcary $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################

# Determine compiler version
SciPrintString("")
include(${SCIMAKE_DIR}/SciCxxFindVersion.cmake)
if (CXX_VERSION)
  SciPrintVar(CXX_VERSION)
else ()
  message(FATAL_ERROR "Could not determine compiler version.")
endif ()

# Set the lib subdir from the Compiler ID and version
if (DEBUG_CMAKE)
  SciPrintVar(CMAKE_CXX_COMPILER_ID)
endif ()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL GNU OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL Clang)
  if (NOT USING_MINGW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -pipe")
  endif ()
  string(SUBSTRING ${CXX_VERSION} 0 1 CXX_MAJOR_VERSION)
  set(CXX_COMP_LIB_SUBDIR gcc${CXX_MAJOR_VERSION})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL Cray)
  string(REGEX REPLACE "\\.[0-9]+-.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  set(CXX_COMP_LIB_SUBDIR cray${CXX_MAJOR_VERSION})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL Intel)
  string(REGEX REPLACE "\\.[0-9]+.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  set(CXX_COMP_LIB_SUBDIR icpc${CXX_MAJOR_VERSION})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL PathScale)
  string(SUBSTRING ${CXX_VERSION} 0 1 CXX_MAJOR_VERSION)
  set(CXX_COMP_LIB_SUBDIR path${CXX_MAJOR_VERSION})
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL PGI)
  string(REGEX REPLACE "\\.[0-9]+-.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  set(CXX_COMP_LIB_SUBDIR pgi${CXX_MAJOR_VERSION})
# Don't automatically include standard library headers.
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --no_using_std")
# Compiler optimization flags set based on "ultra" optimization in
# flags.m4.  Overrides scimake default, since that had -Mipa=fast (no inline).
  set(CMAKE_CXX_FLAGS_RELEASE
    "-fast -O3 -DNDEBUG -Munroll -Minline=levels:5 -Mipa=fast,inline -Mmovnt")
# For a fully-optimized build, set IPA options for linker too
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE
    "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -Mipa=fast,inline")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL XL)
# This should be the basename of the compiler
  string(REGEX REPLACE "\\.[0-9]+.*$" "" CXX_MAJOR_VERSION ${CXX_VERSION})
  string(REGEX REPLACE "^0+" "" CXX_MAJOR_VERSION ${CXX_MAJOR_VERSION})
  get_filename_component(REL_CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} NAME)
# Since we install ben builds in a completely different directory, can
# use same name for CXX_COMP_LIB_SUBDIR
  if (${REL_CMAKE_CXX_COMPILER} MATCHES ".*_r$")
    set(CXX_COMP_LIB_SUBDIR xlC_r${CXX_MAJOR_VERSION})
  else ()
    set(CXX_COMP_LIB_SUBDIR xlC${CXX_MAJOR_VERSION})
  endif ()
  set(SEPARATE_INSTANTIATIONS 1 CACHE BOOL "Whether to separate instantiations -- for correct compilation on xl")
endif ()
SciPrintVar(CXX_COMP_LIB_SUBDIR)

# Look for includes
include(CheckCXXSourceCompiles)
include(CheckIncludeFileCXX)
check_include_file_cxx(sstream HAVE_SSTREAM)
check_include_file_cxx(iostream HAVE_IOSTREAM)

# See whether generally declared statics work
try_compile(HAVE_GENERALLY_DECLARED_STATICS ${PROJECT_BINARY_DIR}/scimake
  ${SCIMAKE_DIR}/trycompile/gendeclstatics.cxx)
set(HAVE_GENERALLY_DECLARED_STATICS ${HAVE_GENERALLY_DECLARED_STATICS} CACHE BOOL "Whether the C++ compiler allows generally declared templated static variables")
if (HAVE_GENERALLY_DECLARED_STATICS)
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/gendeclstatics.cxx compiled.")
  endif ()
else ()
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/gendeclstatics.cxx did not compile.")
  endif ()
endif ()

# See whether std::abs<double> known.
try_compile(HAVE_STD_ABS_DOUBLE ${PROJECT_BINARY_DIR}/scimake
  ${SCIMAKE_DIR}/trycompile/stdabsdbl.cxx)
set(HAVE_STD_ABS_DOUBLE ${HAVE_STD_ABS_DOUBLE} CACHE BOOL "Whether the C++ compiler understands std::abs with double arg")
if (HAVE_STD_ABS_DOUBLE)
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/stdabsdbl.cxx compiled.")
  endif ()
else ()
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/stdabsdbl.cxx did not compile.")
  endif ()
  set(NOT_HAVE_STD_ABS_DOUBLE 1 CACHE BOOL "Define when the C++ compiler does not understand std::abs with double arg")
endif ()

# See whether compiler RTTI typeid is working properly
try_run(RTTI_RUN_RESULT RTTI_COMPILES ${PROJECT_BINARY_DIR}/scimake
  ${SCIMAKE_DIR}/trycompile/checkCompilerRTTI.cxx)
message(STATUS "RTTI_RUN_RESULT = ${RTTI_RUN_RESULT}.")
message(STATUS "RTTI_COMPILES = ${RTTI_COMPILES}.")
set(RTTI_RUN_RESULT ${RTTI_RUN_RESULT} CACHE BOOL "Whether the C++ compiler builds executables that understand run-time type identification.")
set(RTTI_COMPILES ${RTTI_COMPILES} CACHE BOOL "Whether the C++ compiler compiles source using run-time type identification.")
if (RTTI_COMPILES)
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/checkCompilerRTTI.cxx compiled.")
  endif ()
  if (RTTI_RUN_RESULT EQUAL 0)
    set(COMPILER_TYPEID_IS_VALID 1)
    if (DEBUG_CMAKE)
      message(STATUS "Compiler RTTI typeid test passed.")
    endif ()
  elseif ()
    if (DEBUG_CMAKE)
      message(WARNING "Compiler RTTI typeid test did not pass.")
    endif ()
  endif ()
else ()
  if (DEBUG_CMAKE)
    message("${SCIMAKE_DIR}/trycompile/checkCompilerRTTI.cxx did not compile.")
  endif ()
endif ()

include(CheckCXXSourceCompiles)

# Check for iterator being same as pointer
check_cxx_source_compiles(
"
#include <vector>
void f(int* i){}
void f(std::vector<int>::iterator i){}
int main(int argc, char** argv) {return 0;}
"
VECTOR_ITERATOR_IS_NOT_POINTER
)
if (VECTOR_ITERATOR_IS_NOT_POINTER)
  if (DEBUG_CMAKE)
    message(STATUS "std::vector<int>::iterator and int* are not the same.")
  endif ()
else ()
  if (DEBUG_CMAKE)
    message(STATUS "std::vector<int>::iterator and int* are the same.")
  endif ()
  set(VECTOR_ITERATOR_IS_NOT_POINTER 1 CACHE BOOL "Whether std::vector<int>::iterator is the same as int*")
endif ()

# Add in full flags
set(CMAKE_CXX_FLAGS_FULL "${CMAKE_C_FLAGS_FULL}")

# Remove /MD etc for static builds on Windows
if (WIN32 AND NOT MINGW AND NOT BUILD_WITH_SHARED_RUNTIME)
  foreach (flag_var CMAKE_CXX_FLAGS_FULL CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_DEBUG)
    if (NOT (BUILD_SHARED_LIBS OR BUILD_WITH_SHARED_RUNTIME))
      string(REPLACE "/MDd" "" ${flag_var} "${${flag_var}}")
      string(REPLACE "/MD" "" ${flag_var} "${${flag_var}}")
    endif ()
    set(${flag_var} "${${flag_var}} /bigobj")
  endforeach (flag_var)
endif ()

# Check flags
foreach (bld FULL RELEASE RELWITHDEBINFO MINSIZEREL DEBUG)
  SciPrintVar(CMAKE_CXX_FLAGS_${bld})
endforeach ()
SciPrintVar(CMAKE_CXX_FLAGS)

