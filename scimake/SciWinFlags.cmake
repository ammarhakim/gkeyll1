#
# Include for common Windows flags and settings.
#
# $Id: SciWinFlags.cmake 1331 2012-04-24 18:52:42Z mkoch $
#
######################################################################

if (WIN32)
# ICL needs to be defined here(Intel compiler and Visual Studio) for
# trilinos: ml_utils.h
  # add_definitions(-DWIN32 -DICL)
# Windows does not mean ICL!
  add_definitions(-DWIN32)
  set(_USE_MATH_DEFINES 1
    CACHE STRING "Define whether to use math defines(for Windows)")
  if (DEBUG_CMAKE)
    message("CMAKE_C_COMPILER = ${CMAKE_C_COMPILER}.")
  endif ()
  string(REGEX MATCH "^.*icl\\.*" USING_ICL "${CMAKE_C_COMPILER}")
  # MESSAGE("USING_ICL = ${USING_ICL}.")
  string(REGEX MATCH "^.*cl\\.*" USING_CL "${CMAKE_C_COMPILER}")
# Below needs to be fixed to use output of $CC -v
  string(REGEX MATCH "^.*mingw.*" USING_MINGW "${CMAKE_C_COMPILER}")
  # MESSAGE("USING_CL = ${USING_CL}.")
  if (USING_ICL)
    add_definitions(-DICL)
    set(_TIMEVAL_DEFINED 1
        CACHE STRING "Define whether system has timeval(for Windows)")
    # MESSAGE("USING_ICL is true.")
# JRC, 23apr2012: Commenting out this code, as it prevents command-line
# override, and I believe these are defined by scimake.
    if (NOT BUILD_SHARED_LIBS AND NOT BUILD_WITH_SHARED_RUNTIME)
      set(CMAKE_CXX_FLAGS_DEBUG "/MTd /Z7 /Od")
      set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2")
      set(CMAKE_CXX_FLAGS_MINSIZEREL "/MT /O2")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MTd /Z7 /Od")
      set(CMAKE_C_FLAGS_DEBUG "/MTd /Z7 /Od")
      set(CMAKE_C_FLAGS_RELEASE "/MT /O2")
      set(CMAKE_C_FLAGS_MINSIZEREL "/MT /O2")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO "/MTd /Z7 /Od")
    else ()
      # change nothing, should use cmake defaults
    endif ()
    foreach (i DEBUG RELEASE MINSIZERELEASE REWITHDEBINFO)
      set(CMAKE_C_FLAGS_${i} "${CMAKE_C_FLAGS_${i}} /Qstd:c99")
    endforeach ()
  elseif (USING_CL)
    add_definitions(-DCL)
    set(_TIMEVAL_DEFINED 1
        CACHE STRING "Define whether system has timeval(for Windows)")
    if (NOT BUILD_SHARED_LIBS AND NOT BUILD_WITH_SHARED_RUNTIME)
      set(CMAKE_C_FLAGS_DEBUG "/MT /Zi /Ob0 /Od")
      set(CMAKE_C_FLAGS_MINSIZEREL "/MT /O1 /Ob1 /D NDEBUG")
      set(CMAKE_C_FLAGS_RELEASE "/MT /O2 /Ob2 /D NDEBUG")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO "/MT /Zi /O2 /Ob1 /D NDEBUG")
      set(CMAKE_CXX_FLAGS_DEBUG "/MT /Zi /Ob0 /Od")
      set(CMAKE_CXX_FLAGS_MINSIZEREL "/MT /O1 /Ob1 /D NDEBUG")
      set(CMAKE_CXX_FLAGS_RELEASE "/MT /O2 /Ob2 /D NDEBUG")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MT /Zi /O2 /Ob1 /D NDEBUG")
    else ()
      # change nothing, should use cmake defaults
    endif ()
  endif ()
endif ()

