##############################################################################
# 
# SciSquishTestScript.cmake: Commands to be executed by Squish ctests
#
# $Id: SciSquishTestScript.cmake 1363 2012-06-08 18:39:04Z jdelamere $
#
# This script launches a GUI test using Squish.  You should not call
# the script directly; instead, you should access it via the
# macro that is defined in SciSquishMacros.cmake
#
##############################################################################

set(ENV{SQUISH_LIBQTDIR} ${Squish_LIBQTDIR})
set(ENV{SQUISH_USER_SETTINGS_DIR} ${Squish_USER_SETTINGS_DIR})
set(ENV{SQUISH_LICENSEKEY_DIR} ${Squish_LICENSEDIR})

file(WRITE envvars "SQUISH_LIBQTDIR=${Squish_LIBQTDIR}\n")
file(APPEND envvars "SQUISH_USER_SETTINGS_DIR=${Squish_USER_SETTINGS_DIR}\n")
file(APPEND envvars "SQUISH_LICENSEKEY_DIR=${Squish_LICENSEDIR}\n")

if (WIN32)
  set(Bash_START "C:/cygwin/bin/bash.exe")
  set(Bash_ARG "--login")
  set(txcorpdir "$ENV{USERPROFILE}/Documents/txcorp")
  if(IS_DIRECTORY ${scicorpdir})
    file(REMOVE_RECURSE ${scicorpdir})
  endif ()
else ()
  set(Bash_START "")
  set(Bash_ARG "")
  set(txcorpdir "$ENV{HOME}/Documents/txcorp")
  if(IS_DIRECTORY ${scicorpdir})
    file(REMOVE_RECURSE ${scicorpdir})
  endif ()
  if (APPLE)
    set(ENV{DYLD_LIBRARY_PATH} "${Squish_LDLIBPATH}")
    file(APPEND envvars "DYLD_LIBRARY_PATH=${Squish_LDLIBPATH}")    
  else ()
    set(ENV{LD_LIBRARY_PATH} "${Squish_LDLIBPATH}")
    file(APPEND envvars "LD_LIBRARY_PATH=${Squish_LDLIBPATH}")
  endif ()
endif ()


execute_process(
  COMMAND ${Bash_START} ${Bash_ARG} ${Squish_SCRIPT} ${Squish_SERVER_EXECUTABLE} ${Squish_CLIENT_EXECUTABLE} ${Squish_AUT} ${Squish_AUTPATH} ${Ctests_TESTSUITE} ${Ctests_TEST} ${Ctests_LOGFILE}   
  RESULT_VARIABLE ERRTOT
)

# check for an error with running the test
if(NOT "${ERRTOT}" STREQUAL "0")
  message(FATAL_ERROR "Squish Errors: ${ERRTOT}")
endif(NOT "${ERRTOT}" STREQUAL "0")



