######################################################################
#
# SciFindPackage: find includes and libraries of a package
#
# $Id: SciFindPackage.cmake 1285 2012-03-05 20:50:05Z cary $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################
#
# Variables which can be used to change the behavior of this script
#
# Where "scipkgreg" is the regularized package name
#(package name with all "." and "-" replaced with "_")
# Where "scipkguc" is the UPPERCASE regularized package name
#
#  DEBUG_CMAKE - if true, outputs verbose debugging information
#(default false)
#
#  ENABLE_${scipkguc} - if false, will not search for package
#(default true)
#  DISABLE_${scipkguc} - if true, sets ENABLE_${scipkguc} to false
#
#  ${scipkgreg}_FIND_QUIETLY - if true, will succeed silently
#(default - not defined, which is treated as false)
#
#  ${scipkgreg}_FIND_REQUIRED - if true, will issue a fatal error if
#    package not found
#(default - not defined, which is treated as false)
#
#  ${scipkguc}_ROOT_DIR, ${scipkguc}_DIR - search directory hints
#
#  SUPRA_SEARCH_PATH - used to specify various top-level search directories
#
######################################################################
#
#  Variables which will be defined by this script
#  In the event that a variable has no valid value, it will be set to
#   "${varname}-NOTFOUND"
#
#  NOTE: Some variables use scipkguc, and some use scipkgreg.
#    Caveat emptor.
#
#    ${scipkguc}_FOUND - true if package found
#
#  EXECUTABLES:
#    ${scipkgreg}_EXECUTABLES - list of found executables, including full
#      path to each.  Only defined for specifically requested executables.
#    ${scipkgreg}_yyy - full name & path to executable "yyy"
#      Only defined for specifically requested executables.
#
#  FILES:
#    ${scipkgreg}_FILES - list of found files, including a full path to each
#    ${scipkgreg}_yyy - full name and path to file "yyy"
#
#  MODULES:
#    ${scipkgreg}_MODULE_DIRS - a list of all module directories found
#    ${scipkgreg}_yyy_MOD - full path to individual module yyy.{mod,MOD}
#      Only defined for specifically requested modules
#
#  HEADERS:
#    ${scipkgreg}_INCLUDE_DIRS - a list of all include directories found
#    ${scipkgreg}_yyy - full path to individual header yyy
#      Only defined for specifically requested headers
#
#  LIBRARIES:
#    ${scipkgreg}_{library_name}_LIBRARY - full path to the individual
#      $library_name library. Only defined for specifically requested libraries.
#    ${scipkgreg}_LIBRARIES - list of found libraries, including full path to
#      each.  Only defined for specifically requested libraries.
#   ${scipkgreg}_STLIBS - list of all found static libraries.
#      Only defined for libraries existing in ${scipkgreg}_LIBRARIES
#   ${scipkgreg}_DLLS - windows only, list of all found dlls
#      Only defined for libraries existing in ${scipkgreg}_LIBRARIES
#
#######################################################################

# SciGetStaticLibs
#
# Given a list of libraries, create a new list, where for any
# dynamic library with a static library next to it, the new list
# contains the static library.  Otherwise it contains the original
# dynamic library.
#
# Args:
# origlibs: the original libraries
# statlibsvar: the variable holding the static libs when found
#
FUNCTION(SciGetStaticLibs origlibs statlibsvar)
  # set(origlibs ${${origlibsvar}})
  set(statlibs)
  if (DEBUG_CMAKE)
    message(STATUS "[SciFindPackage]: original libs = ${origlibs}.")
    message(STATUS "[SciFindPackage]: statlibsvar = ${statlibsvar}.")
  endif ()
  foreach (lib ${origlibs})
    set(newlib ${lib})
    set(havestlib FALSE)
    if (${lib} MATCHES "\\.a$" OR ${lib} MATCHES "\\.lib$")
      if (DEBUG_CMAKE)
        message(STATUS "${lib} is a static library.")
      endif ()
      set(havestlib TRUE)
    endif ()
# If not static, try replacing suffix
    if (NOT ${havestlib})
      foreach (sfx so;dylib;dll)
        if (${lib} MATCHES "\\.${sfx}$")
          if (DEBUG_CMAKE)
            message(STATUS "${lib} is a shared library.")
          endif ()
          get_filename_component(libdir ${lib}/.. REALPATH)
# NAME_WE takes from first .  We generally need from last
          get_filename_component(libname "${lib}" NAME)
          string(REGEX REPLACE "\\.[^\\.]*$" "" libname "${libname}")
          if (${sfx} STREQUAL "dll")
            set(newsfx "lib")
          else ()
            set(newsfx "a")
          endif ()
          if (EXISTS ${libdir}/${libname}.${newsfx})
            set(newlib ${libdir}/${libname}.${newsfx})
            set(havestlib TRUE)
          endif ()
          break()
        endif ()
      endforeach ()
    endif ()
# If still do not have static, try pulling out of library flags
    if (${havestlib})
      list(APPEND statlibs ${newlib})
    elseif (${lib} MATCHES "^-L")
      if (DEBUG_CMAKE)
        message(STATUS "${lib} defined by flags.")
      endif ()
      # set(libargs ${lib}) # Create a list
      string(REPLACE " " ";" libargs "${lib}")
      if (DEBUG_CMAKE)
        message(STATUS "libargs = ${libargs}.")
      endif ()
      set(libdirargs)
      set(libdirs)
      set(libnames)
      foreach (libarg ${libargs})
        if (DEBUG_CMAKE)
          message(STATUS "libarg = ${libarg}.")
        endif ()
        if (${libarg} MATCHES "^-L")
          list(APPEND libdirargs ${libarg})
          string(REGEX REPLACE "^-L" "" libdir ${libarg})
          list(APPEND libdirs ${libdir})
        elseif (${libarg} MATCHES "^-l")
          string(REGEX REPLACE "^-l" "" libname ${libarg})
          list(APPEND libnames ${libname})
        endif ()
      endforeach ()
      if (DEBUG_CMAKE)
        message(STATUS "libdirs = ${libdirs}.")
        message(STATUS "libnames = ${libnames}.")
      endif ()
      set(stlibs)
      foreach (libname ${libnames})
        set(lib)
        foreach (libdir ${libdirs})
          if (EXISTS ${libdir}/lib${libname}.a)
            set(lib ${libdir}/lib${libname}.a)
            break()
          elseif (EXISTS ${libdir}/${libname}.lib)
            set(lib ${libdir}/${libname}.lib)
            break()
          endif ()
        endforeach ()
        if (lib)
          list(APPEND stlibs ${lib})
        else ()
          list(APPEND stlibs "${libdirargs} -l${libname}")
        endif ()
      endforeach ()
      list(APPEND statlibs ${stlibs})
    else ()
# If still not static, look in system dirs
      foreach (libdir /usr/lib64;/usr/lib)
        if (EXISTS ${libdir}/lib${lib}.a)
          set(newlib ${libdir}/lib${lib}.a)
          set(havestlib TRUE)
          break()
        endif ()
      endforeach ()
      list(APPEND statlibs ${newlib})
    endif ()
  endforeach ()
  if (DEBUG_CMAKE)
    message(STATUS "[SciFindPackage]: static libs = ${statlibs}.")
  endif ()
  set(${statlibsvar} ${statlibs} PARENT_SCOPE)
ENDFUNCTION()

# SciGetRealDir
#
# Given a directory name, find the real directory name after
#  first trying to resolve it as a shortcut if on Windows,
#  then, if that does not work, trying to resolve it as a soft
#  link.
#
# Args:
#  canddir: the variable holding the directory candidate
#  realdirvar: the variable holding the directory after following any shortcuts
#
FUNCTION(SciGetRealDir canddir realdirvar)
  # MESSAGE("realdirvar = ${realdirvar}.")
  set(${realdirvar})
  if (canddir)
    if (DEBUG_CMAKE)
      message(STATUS "Finding the directory for ${canddir}.")
    endif ()
# Follow an existing Windows shortcut
    if (WIN32)
      if (EXISTS ${canddir}.lnk)
        exec_program(readshortcut ARGS ${canddir}.lnk
          OUTPUT_VARIABLE idir
          return_VALUE ires)
        if (ires)
          message(FATAL_ERROR "readshortcut error on ${canddir}.lnk.")
        endif ()
        exec_program(cygpath ARGS -am ${idir} OUTPUT_VARIABLE rd)
      else ()
        set(rd ${canddir})
      endif ()
    else ()
       get_filename_component(rd ${canddir} REALPATH)
    endif ()
  endif ()
  set(${realdirvar} ${rd} PARENT_SCOPE)
ENDFUNCTION()

#
# SciFindPackage
#
# Args:
# PACKAGE: the prefix for the defined variables
# INSTALL_DIRS: the names for the installation subdirectory.  Defaults
#   to lower cased scipkgname
# EXECUTABLES: the executables to look for
# HEADERS: the header files to look for.
# LIBRARIES: the libraries
# MODULES: Fortran module files
# EXECUTABLE_SUBDIRS: executable subdirs
# INCLUDE_SUBDIRS: include subdirectories
# LIBRARY_SUBDIRS: library subdirectories
#
# NOTE: One cannot reset calling variables
# NOTE: lists should be delimited by semicolons.
#(which is the default format used by scimake)
#
include(CMakeParseArguments)
macro(SciFindPackage)
  CMAKE_PARSE_ARGUMENTS(TFP
    ""
    "PACKAGE;INSTALL_DIR"
    "INSTALL_DIRS;EXECUTABLES;HEADERS;LIBRARIES;FILES;MODULES;EXECUTABLE_SUBDIRS;INCLUDE_SUBDIRS;MODULE_SUBDIRS;LIBRARY_SUBDIRS;FILE_SUBDIRS"
    ${ARGN}
  )

  # This message is purposefully NOT a STATUS message
  # To provide more readable output
  message("")
  message("--------- SciFindPackage looking for ${TFP_PACKAGE} ---------")

  if (DEBUG_CMAKE)
    message(STATUS "Using verbose options")
    message(STATUS "SciFindPackage called with arguments:
    PACKAGE = ${TFP_PACKAGE}
    INSTALL_DIR = ${TFP_INSTALL_DIR}
    INSTALL_DIRS = ${TFP_INSTALL_DIRS}
    EXECUTABLES = ${TFP_EXECUTABLES}
    HEADERS = ${TFP_HEADERS}
    LIBRARIES = ${TFP_LIBRARIES}
    MODULES = ${TFP_MODULES}
    EXECUTABLE_SUBDIRS = ${TFP_EXECUTABLE_SUBDIRS}
    INCLUDE_SUBDIRS = ${TFP_INCLUDE_SUBDIRS}
    LIBRARY_SUBDIRS = ${TFP_LIBRARY_SUBDIRS}")
  endif ()

# Construct various names(upper/lower case) for package
  string(REGEX REPLACE "[./-]" "_" scipkgreg ${TFP_PACKAGE})
# scipkgreg is the regularized package name
  string(TOUPPER ${scipkgreg} scipkguc)
  if (DEBUG_CMAKE)
    message(STATUS "ENABLE_${scipkguc} = ${ENABLE_${scipkguc}}")
  endif ()
  if (DISABLE_${scipkguc})
    set(ENABLE_${scipkguc} FALSE)
  endif ()
  if (NOT DEFINED ENABLE_${scipkguc})
    set(ENABLE_${scipkguc} TRUE)
  endif ()
  if (NOT ENABLE_${scipkguc})
    message(STATUS "Disabling ${scipkgreg}.")
    set(${scipkguc}_FOUND FALSE)
    if (DEBUG_CMAKE)
      message(STATUS "${scipkguc}_FOUND set to FALSE.")
    endif ()
    return()
  endif ()
  string(TOLOWER ${scipkgreg} scipkglc)
  set(scipkginst ${TFP_INSTALL_DIR} ${TFP_INSTALL_DIRS})
  if (NOT scipkginst)
    set(scipkginst ${scipkginst} ${scipkglc})
  endif ()
  if (DEBUG_CMAKE)
    message(STATUS "scipkginst = ${scipkginst}.")
  endif ()

# Find the possible directories and put them into path
# Order: command-line define, environment variable, supra-search-path subdirs,
# supra-search-path dirs
  if (DEBUG_CMAKE)
    message("SciFindPackage] SUPRA_SEARCH_PATH = ${SUPRA_SEARCH_PATH}")
  endif ()
  set(scipath)
# Command-line define overrides all.
  if (${scipkg}_ROOT_DIR)
    SciGetRealDir(${${scipkg}_ROOT_DIR} scipath)
  elseif (${scipkg}_DIR)
# JRC 20120617: Remove this July 31, 2012.
    message(WARNING "Use of ${scipkg}_DIR define is deprecated.  Please use ${scipkg}_ROOT_DIR")
    SciGetRealDir(${${scipkg}_DIR} scipath)
  endif ()
  if (NOT DEFINED ${scipath})
    if ($ENV{${scipkg}_ROOT_DIR})
      SciGetRealDir($ENV{${scipkg}_ROOT_DIR} scipath)
    elseif (${scipkg}_DIR)
# JRC 20120617: Remove this July 31, 2012.
      message(WARNING "Use of ${scipkg}_DIR environment variable is deprecated.  Please use ${scipkg}_ROOT_DIR")
      SciGetRealDir($ENV{${scipkg}_DIR} scipath)
    endif ()
  endif ()
# Next try environment variable
# Supra-search-path dirs
  foreach (instdir ${scipkginst})
    foreach (spdir ${SUPRA_SEARCH_PATH})
      set(idir ${spdir}/${instdir})
      if (EXISTS ${idir})
        SciGetRealDir(${idir} scidir)
        set(scipath ${scipath} ${scidir})
      endif ()
    endforeach (spdir ${SUPRA_SEARCH_PATH})
  endforeach ()
# Supra-search-path dirs
  foreach (spdir ${SUPRA_SEARCH_PATH})
    set(idir ${spdir})
    if (EXISTS ${idir})
      SciGetRealDir(${idir} scidir)
      set(scipath ${scipath} ${scidir})
    endif ()
  endforeach (spdir ${SUPRA_SEARCH_PATH})
# Any found?
  list(LENGTH scipath scipathlen)
  if (DEBUG_CMAKE)
    if (NOT scipathlen)
      message(FATAL_ERROR "scipath is empty.")
    else ()
      message(STATUS "scipath = ${scipath}")
    endif ()
  endif ()

#######################################################################
#
# Look for EXECUTABLES
# Variables defined:
#   XXX_YYY - CACHED
#     Where to find the YYY tool that comes with XXX.
#   XXX_EXECUTABLES - CACHED
#     List of all executables found for package XXX.
#
#######################################################################

# Create the search paths
  string(LENGTH "${TFP_EXECUTABLES}" sciexecslen)
  if (${sciexecslen})
    string(LENGTH "${TFP_EXECUTABLE_SUBDIRS}" scilen)
    if (${scilen})
      set(sciexecsubdirs ${TFP_EXECUTABLE_SUBDIRS})
    else ()
# Here we use the default search directory "bin"
      set(sciexecsubdirs bin)
    endif ()

    if (DEBUG_CMAKE)
      message(STATUS "Looking for executables under ${scipath} with subdirs, ${sciexecsubdirs}")
    endif ()

# Clear the list of executable dirs
# As we loop over the executables, we will add every new directory to the list
# And afterwards, we'll remove duplicates
    set(${scipkgreg}_EXECUTABLES)
    set(${scipkgreg}_FOUND_SOME_EXECUTABLE FALSE)

# Loop over list of executables and try to find each one
    foreach (sciexec ${TFP_EXECUTABLES})

# Create the variable that's specific to this executable
      string(REGEX REPLACE "[./-]" "_" sciexecvar ${sciexec})
      set(sciexecvar ${scipkgreg}_${sciexecvar})
      if (DEBUG_CMAKE)
        message(STATUS "Calling FIND_PROGRAM with sciexec = ${sciexec}")
        message(STATUS "Calling FIND_PROGRAM with paths = ${scipath}")
        message(STATUS "Calling FIND_PROGRAM with path_suffixes ${sciexecsubdirs}")
      endif ()

# First look in specified path
      find_program(${sciexecvar}
        "${sciexec}"
        PATHS ${scipath}
        PATH_SUFFIXES "${sciexecsubdirs}"
        NO_DEFAULT_PATH
        DOC "Path to the ${sciexec} executable"
      )

# If not found, try again with default paths
      if (NOT ${sciexecvar})
        if (DEBUG_CMAKE)
          message(STATUS "Failed to find ${sciexec} in search path, trying default paths.")
        endif ()

        find_program(${sciexecvar}
          ${sciexec}
          PATHS ${scipath}
          PATH_SUFFIXES "${sciexecsubdirs}"
          DOC "Path to the ${sciexec} executable"
        )
      endif ()

      if (DEBUG_CMAKE)
        message(STATUS "${sciexecvar} is ${${sciexecvar}}")
      endif ()

      if (${${sciexecvar}} MATCHES NOTFOUND)
# The WARNING option will actually give a scimake stack trace
# And I didn't want that
# So instead, use NO option, and start the string with WARNING
        message("WARNING - Unable to locate executable ${sciexec}")
      else ()
        set(${scipkgreg}_FOUND_SOME_EXECUTABLE TRUE)
# Add to list of all executables for this package
        set(${scipkgreg}_EXECUTABLES
          ${${scipkgreg}_EXECUTABLES}
          ${${sciexecvar}}
        )
      endif ()
    endforeach (sciexec ${TFP_EXECUTABLES})

# Clean up the list of all executable directories
    if (${scipkgreg}_FOUND_SOME_EXECUTABLE)
      list(REMOVE_DUPLICATES "${scipkgreg}_EXECUTABLES")
    endif ()

# Print results if in debug mode
    if (DEBUG_CMAKE)
      message(STATUS "List of all executables for ${scipkgreg}:")
      message(STATUS "${scipkgreg}_EXECUTABLES = ${${scipkgreg}_EXECUTABLES}")
    endif ()

# Cache the completed list
    set(${scipkgreg}_EXECUTABLES
      ${${scipkgreg}_EXECUTABLES}
      CACHE STRING "List of all executables for ${scipkgreg}"
    )
  endif ()

#######################################################################
#
# Look for FILES.
# Like EXECUTABLES, but using find_file, instead of find_program
# Variables defined:
#   XXX_YYY - CACHED
#     Where to find the YYY tool that comes with XXX.
#   XXX_FILES - CACHED
#     List of all executables found for package XXX.
#
#######################################################################

# Create the search paths
  string(LENGTH "${TFP_FILES}" scifileslen)
  if (${scifileslen})
    string(LENGTH "${TFP_FILE_SUBDIRS}" scilen)
    if (${scilen})
      set(scifilesubdirs ${TFP_FILE_SUBDIRS})
    else ()
# Here we use the default search directory "bin"
      # set(scifilesubdirs .)
    endif ()

    if (DEBUG_CMAKE)
      message(STATUS "Looking for files under ${scipath} with subdirs, ${scifilesubdirs}")
    endif ()

# Clear the list of file dirs
# As we loop over the files, we will add every new directory to the list
# And afterwards, we'll remove duplicates
    set(${scipkgreg}_FILES)
    set(${scipkgreg}_FOUND_SOME_FILE FALSE)

# Loop over list of files and try to find each one
    foreach (scifile ${TFP_FILES})

# Create the variable that's specific to this file
      string(REGEX REPLACE "[./-]" "_" scifilevar ${scifile})
      set(scifilevar ${scipkgreg}_${scifilevar})
      if (DEBUG_CMAKE)
        message(STATUS "Calling FIND_FILE with scifile = ${scifile}")
        message(STATUS "Calling FIND_FILE with paths = ${scipath}")
        message(STATUS "Calling FIND_FILE with path_suffixes ${scifilesubdirs}")
      endif ()

# First look in specified path
      find_file(${scifilevar}
        "${scifile}"
        PATHS ${scipath}
        PATH_SUFFIXES ${scifilesubdirs}
        NO_DEFAULT_PATH
        DOC "Path to the ${scifile} file"
      )

# If not found, try again with default paths
      if (NOT ${scifilevar})
        if (DEBUG_CMAKE)
          message(STATUS "Failed to find ${scifile} in search path, trying default paths.")
        endif ()

        find_file(${scifilevar}
          ${scifile}
          PATHS ${scipath}
          PATH_SUFFIXES ${scifilesubdirs}
          DOC "Path to the ${scifile} file"
        )
      endif ()

      if (DEBUG_CMAKE)
        message(STATUS "${scifilevar} is ${${scifilevar}}")
      endif ()

      if (${${scifilevar}} MATCHES NOTFOUND)
# The WARNING option will actually give a scimake stack trace
# And I didn't want that
# So instead, use NO option, and start the string with WARNING
        message("WARNING - Unable to locate file ${scifile}")
      else ()
        set(${scipkgreg}_FOUND_SOME_FILE TRUE)
# Add to list of all files for this package
        set(${scipkgreg}_FILES
          ${${scipkgreg}_FILES}
          ${${scifilevar}}
        )
      endif ()
    endforeach (scifile ${TFP_FILES})

# Clean up the list of all file directories
    if (${scipkgreg}_FOUND_SOME_FILE)
      list(REMOVE_DUPLICATES "${scipkgreg}_FILES")
    endif ()

# Print results if in debug mode
    if (DEBUG_CMAKE)
      message(STATUS "List of all files for ${scipkgreg}:")
      message(STATUS "${scipkgreg}_FILES = ${${scipkgreg}_FILES}")
    endif ()

# End of this block
  endif ()

#######################################################################
#
# Look for MODULES
# Like FILES, but append module suffix
# Variables defined:
#   XXX_YYY - CACHED
#     Where to find the YYY tool that comes with XXX.
#   XXX_MODULES - CACHED
#     List of all executables found for package XXX.
#
#######################################################################

# Create the search paths
  string(LENGTH "${TFP_MODULES}" scimoduleslen)
  if (${scimoduleslen} AND NOT SCI_FC_MODULE_SUFFIX)
    message(WARNING "[SciFindPackage.cmake] Cannot find Fortran module files, as SCI_FC_MODULE_SUFFIX is not defined with a call to SciFortranChecks.cmake.")
  elseif (${scimoduleslen})
    string(LENGTH "${TFP_MODULE_SUBDIRS}" scilen)
    if (${scilen})
      set(scimodulesubdirs ${TFP_MODULE_SUBDIRS})
    else ()
# Default subdirectory
      set(scimodulesubdirs include lib)
    endif ()

    if (DEBUG_CMAKE)
      message(STATUS "Looking for modules under ${scipath} with subdirs, ${scimodulesubdirs}")
    endif ()

# Clear the list of file dirs
# As we loop over the files, we will add every new directory to the list
# And afterwards, we'll remove duplicates
    set(${scipkgreg}_MODULES)
    set(${scipkgreg}_MODULE_DIRS)
    set(${scipkgreg}_FOUND_SOME_MODULE FALSE)

# Loop over list of files and try to find each one
    foreach (scimodule ${TFP_MODULES})

# Capitalize as needed
      if (SCI_FC_MODULENAME_CAPITALIZED)
        string(TOUPPER ${scimodule} scimodfile)
      else ()
        set(scimodfile ${scimodule})
      endif ()
      set(scimodfile "${scimodfile}.${SCI_FC_MODULE_SUFFIX}")

# Create the variable that's specific to this file
      string(REGEX REPLACE "[./-]" "_" scimodulevar ${scimodule})
      set(scimodulevar ${scipkgreg}_${scimodulevar}_MOD)
      if (DEBUG_CMAKE)
        message(STATUS "Calling FIND_FILE with scimodfile = ${scimodfile}")
        message(STATUS "Calling FIND_FILE with paths = ${scipath}")
        message(STATUS "Calling FIND_FILE with path_suffixes ${scimodulesubdirs}")
      endif ()

# First look in specified path
      find_file(${scimodulevar}
        "${scimodfile}"
        PATHS ${scipath}
        PATH_SUFFIXES ${scimodulesubdirs}
        NO_DEFAULT_PATH
        DOC "Path to the ${scimodfile} file"
      )

# If not found, try again with default paths
      if (NOT ${scimodulevar})
        if (DEBUG_CMAKE)
          message(STATUS "Failed to find ${scimodfile} in search path, trying default paths.")
        endif ()

        find_file(${scimodulevar}
          ${scimodfile}
          PATHS ${scipath}
          PATH_SUFFIXES ${scimodulesubdirs}
          DOC "Path to the ${scimodfile} file"
        )
      endif ()

      if (DEBUG_CMAKE)
        message(STATUS "${scimodulevar} is ${${scimodulevar}}")
      endif ()

      if (${${scimodulevar}} MATCHES NOTFOUND)
# The WARNING option will actually give a scimake stack trace
# And I didn't want that
# So instead, use NO option, and start the string with WARNING
        message("WARNING - Unable to locate file ${scimodule}")
      else ()
        set(${scipkgreg}_FOUND_SOME_MODULE TRUE)
# Add to list of all files for this package
        set(${scipkgreg}_MODULES
          ${${scipkgreg}_MODULES}
          ${${scimodulevar}}
        )
        set(${scimodulevar} ${${scimodulevar}}
          CACHE FILEPATH "${scimodule} file."
        )
        get_filename_component(dir ${${scimodulevar}}/.. REALPATH)
        set(${scipkgreg}_MODULE_DIRS ${${scipkgreg}_MODULE_DIRS} ${dir})
      endif ()
    endforeach (scimodule ${TFP_MODULES})

# Clean up the list of all module directories
    if (${scipkgreg}_FOUND_SOME_MODULE)
      list(REMOVE_DUPLICATES "${scipkgreg}_MODULES")
      list(REMOVE_DUPLICATES "${scipkgreg}_MODULE_DIRS")
    endif ()

# Print results if in debug mode
    if (DEBUG_CMAKE)
      message(STATUS "List of all modules for ${scipkgreg}:")
      message(STATUS "${scipkgreg}_MODULES = ${${scipkgreg}_MODULES}")
    endif ()

# End of this block
  endif ()

##########################################################################
# Look for HEADERS
# Finding none is okay(e.g. zlib has no headers)
# Will set:
#   XXX_INCLUDE_DIRS - NOT CACHED
#     The final set of include directories listed in one variable.
#   XXX_yy_h - NOT CACHED
#     Where to find xxx_yy.h
##########################################################################

  string(LENGTH "${TFP_HEADERS}" scihdrslen)
  if (${scihdrslen})
# Find include subdirs
    string(LENGTH "${TFP_INCLUDE_SUBDIRS}" scilen)
    if (${scilen})
      set(sciincsubdirs ${TFP_INCLUDE_SUBDIRS})
    else ()
      set(sciincsubdirs include)
    endif ()

    if (DEBUG_CMAKE)
      message(STATUS "Looking for headers under ${scipath} with subdirs, ${sciincsubdirs}")
    endif ()

# Clear the variables
    set(${scipkgreg}_INCLUDE_DIRS)
    set(${scipkgreg}_FOUND_SOME_HEADER FALSE)

# Look for each header in turn
    foreach (scihdr ${TFP_HEADERS})

# Create variable name from cleaned name
      string(REGEX REPLACE "[./-]" "_" scihdrvar ${scihdr})
      set(scihdrvar ${scipkgreg}_${scihdrvar})
      set(scihdrdirvar ${scihdrvar}_INCLUDE_DIR)

# First look in specified path
      find_path(${scihdrdirvar}
        ${scihdr}
        PATHS ${scipath}
        PATH_SUFFIXES ${sciincsubdirs}
        NO_DEFAULT_PATH)

# If not found, also look in default paths
      if (NOT ${scihdrdirvar})
        find_path(${scihdrdirvar}
          ${scihdr}
          PATHS ${scipath}
          PATH_SUFFIXES ${sciincsubdirs})
      endif ()

# Found or not?
      if (${scihdrdirvar})
# Create header variable and push into cache
        set(${scihdrdirvar} ${${scihdrdirvar}}
          CACHE FILEPATH "Directory containing ${scihdr}."
        )
        if (DEBUG_CMAKE)
           message(STATUS "${scihdrdirvar} = ${${scihdrdirvar}}.")
        endif ()
        set(${scipkgreg}_FOUND_SOME_HEADER TRUE)
        set(${scihdrvar} ${${scihdrdirvar}}/${scihdr}
          CACHE FILEPATH "${scihdr} file."
        )
# Add the directory of this header to the list of all include directories
        set(${scipkgreg}_INCLUDE_DIRS ${${scipkgreg}_INCLUDE_DIRS}
          ${${scihdrdirvar}}
        )
      else ()
        message("WARNING - Unable to find header: ${scihdr}")
      endif ()

    endforeach (scihdr ${TFP_HEADERS})

# Clean up and save list of include directories into the cache
    list(LENGTH ${scipkgreg}_INCLUDE_DIRS sciinclistlen)
    if (${sciinclistlen})
      list(REMOVE_DUPLICATES ${scipkgreg}_INCLUDE_DIRS)
# scimake conventions are that this is not a cache variable
      # set(${scipkgreg}_INCLUDE_DIRS
        # ${${scipkgreg}_INCLUDE_DIRS}
        # CACHE STRING "All include directories for ${scipkgreg}"
      # )
    else ()
      message(WARNING "None of ${TFP_HEADERS} found.  Define ${scipkguc}_DIR to find them.")
    endif ()
  endif ()

###########################################################################
#
# Look for LIBRARIES
# Finding none is fatal, but finding a subset is okay as allows us
# to look for a maximum set, as needed for Trilinos.
# XXX_LIBRARIES - NOT CACHED
#   The libraries to link against to use XXX. These should include full paths.
# XXX_YY_LIBRARY -
#   Name of YY library that is part of the XXX system.
#
###########################################################################

# Build list of search directories
  string(LENGTH "${TFP_LIBRARIES}" scilibslen)
  if (${scilibslen})
# Add in user-supplied subdirs
    string(LENGTH "${TFP_LIBRARY_SUBDIRS}" scilen)
    if (${scilen})
      set(scilibsubdirs ${TFP_LIBRARY_SUBDIRS})
    else ()
# Default subdirectory
      set(scilibsubdirs lib)
    endif ()

# Does anyone use this?  Should either move to argument list
# or remove entirely.  TODO - marc
# If LIBRARY_DIRS specified, add that to the front of the path
    if (${scipkgreg}_LIBRARY_DIRS)
      set(scilibsubdirs . ${scilibsubdirs})
      set(scipath ${${scipkgreg}_LIBRARY_DIRS} ${scipath})
    endif ()

    if (DEBUG_CMAKE)
      message(STATUS "Looking for libraries under ${scipath} with subdirs, ${scilibsubdirs}")
    endif ()

# Clear variables
    set(${scipkgreg}_LIBRARIES)
    set(${scipkgreg}_LIBRARY_DIRS)
    set(${scipkgreg}_LIBRARY_NAMES)
    set(${scipkgreg}_FOUND_SOME_LIBRARY FALSE)

# Look for each requested library
    foreach (scilib ${TFP_LIBRARIES})
# Build variable name of the form XXX_yyy_LIBRARY
      string(REGEX REPLACE "[./-]" "_" scilibvar ${scilib})
      set(scilibdirvar ${scipkgreg}_${scilibvar}_LIBRARY_DIR)
      unset(${scilibdirvar})
      set(scilibvar ${scipkgreg}_${scilibvar}_LIBRARY)
      unset(${scilibvar})

# First look in defined paths
      find_library(${scilibvar}
        ${scilib}
        PATHS ${scipath}
        PATH_SUFFIXES ${scilibsubdirs}
        NO_DEFAULT_PATH
      )
      if (DEBUG_CMAKE)
        message(STATUS "After initial search, ${scilibvar} = ${${scilibvar}}.")
      endif ()

# If not found, try again in default paths
      if (NOT ${scilibvar})
        find_library(${scilibvar}
          ${scilib}
          PATHS ${scipath}
          PATH_SUFFIXES ${scilibsubdirs}
        )
        if (DEBUG_CMAKE)
          message(STATUS "After second search, ${scilibvar} = ${${scilibvar}}.")
        endif ()
      endif ()

# Add to list of all libraries(if found)
      if (${scilibvar})
        if (DEBUG_CMAKE)
          message(STATUS "Found library ${scilib}: ${scilibvar} = ${${scilibvar}}")
        endif ()
        set(${scipkgreg}_FOUND_SOME_LIBRARY TRUE)
        set(${scipkgreg}_LIBRARIES
            ${${scipkgreg}_LIBRARIES}
            ${${scilibvar}}
        )
        get_filename_component(${scilibdirvar} ${${scilibvar}}/.. REALPATH)
        if (DEBUG_CMAKE)
          message(STATUS "${scilibdirvar} = ${${scilibdirvar}}")
        endif ()
        list(APPEND ${scipkgreg}_LIBRARY_DIRS ${${scilibdirvar}})
# NAME_WE removes from first '.'.  We need from last.
        get_filename_component(libname ${${scilibvar}} NAME)
        string(REGEX REPLACE "\\.[^\\.]*$" "" libname "${libname}")
        string(REGEX REPLACE "^lib" "" libname ${libname})
        list(APPEND ${scipkgreg}_LIBRARY_NAMES ${libname})
      else ()
# Yes, this uses the string WARNING instead of the macro argument WARNING
        message("WARNING - Library ${scilib} not found.")
      endif ()
# Add both to the cache
      # set(${scilibvar}
        # ${${scilibvar}}
        # CACHE FILEPATH "Location of ${scilib}"
      # )
      set(${scilibdirvar}
        ${${scilibdirvar}}
        CACHE PATH "Directory of ${scilib}"
      )
    endforeach (scilib ${TFP_LIBRARIES})
    if (DEBUG_CMAKE)
      message(STATUS "${scipkgreg}_LIBRARY_DIRS = ${${scipkgreg}_LIBRARY_DIRS}")
    endif ()

# Clean up and commit variables to cache
    list(LENGTH ${scipkgreg}_LIBRARIES sciliblistlen)
    if (NOT ${sciliblistlen})
      message("WARNING - None of the libraries, ${TFP_LIBRARIES}, found.  Define ${scipkguc}_DIR to find them.")
    else ()
      list(REMOVE_DUPLICATES "${scipkgreg}_LIBRARIES")
      list(REMOVE_DUPLICATES "${scipkgreg}_LIBRARY_DIRS")
# The first dir is the library dir
      list(GET "${scipkgreg}_LIBRARY_DIRS" 0 ${scipkgreg}_LIBRARY_DIR)
      if (NOT DEFINED ${scipkgreg}_DIR)
        get_filename_component(${scipkgreg}_DIR ${${scipkgreg}_LIBRARY_DIR}/.. REALPATH)
      endif ()
      if (DEBUG_CMAKE)
        message(STATUS "${scipkgreg}_DIR = ${${scipkgreg}_DIR}")
        message(STATUS "${scipkgreg}_LIBRARIES = ${${scipkgreg}_LIBRARIES}")
        message(STATUS "${scipkgreg}_LIBRARY_DIRS = ${${scipkgreg}_LIBRARY_DIRS}")
      endif ()
    endif ()

# Find static libraries
    SciGetStaticLibs("${${scipkgreg}_LIBRARIES}" ${scipkgreg}_STLIBS)

  endif ()

#################################################################
#
# On windows, we want to look for the dlls
# XXX_DLLS - CACHED
#   List of all found dlls including full path to each
# XXX_yyy_DLLS
#   Path to individual yyy library from XXX package
#  DOES NOT EXIST YET - TODO
#
#################################################################

  if (WIN32)
    if (DEBUG_CMAKE)
      message(STATUS "Looking for the DLL counterparts to all found libraries.")
    endif ()

# Clear variables
    set(${scipkgreg}_DLLS)
    set(${scipkgreg}_FOUND_SOME_DLL FALSE)

# Look for a dll for each found library
    foreach (foundlib ${${scipkgreg}_LIBRARIES})
      # get_filename_component(librootname ${foundlib} NAME_WE)
# NAME_WE removes from first '.'.  We need from last.
      get_filename_component(librootname ${foundlib} NAME)
      string(REGEX REPLACE "\\.[^\\.]*$" "" librootname "${librootname}")
# This assumes that dll is in "./lib/../bin"
      get_filename_component(libdir ${foundlib}/.. REALPATH)
      get_filename_component(dlldir1 "${libdir}/../bin" REALPATH)
      get_filename_component(dlldir2 "${libdir}/.." REALPATH)
      if (DEBUG_CMAKE)
        message(STATUS "Looking for DLL counterpart to ${foundlib} in ${dlldir1}/${librootname}.dll")
      endif ()
      if (EXISTS ${dlldir1}/${librootname}.dll)
        if (DEBUG_CMAKE)
          message(STATUS "Found DLL ${dlldir1}/${librootname}.dll")
        endif ()
        set(${scipkgreg}_DLLS ${${scipkgreg}_DLLS} ${dlldir1}/${librootname}.dll)
      else ()
        if (DEBUG_CMAKE)
          message(STATUS "Second chance: looking for DLL in ${dlldir2}/${librootname}.dll")
        endif ()
        if (EXISTS ${dlldir2}/${librootname}.dll)
          if (DEBUG_CMAKE)
            message(STATUS "Found DLL ${dlldir2}/${librootname}.dll")
          endif ()
          set(${scipkgreg}_DLLS ${${scipkgreg}_DLLS} ${dlldir2}/${librootname}.dll)
        else ()
          if (DEBUG_CMAKE)
            message(STATUS "Third chance: looking for DLL in ${libdir}/${librootname}.dll")
          endif ()
          if (EXISTS ${libdir}/${librootname}.dll)
            if (DEBUG_CMAKE)
              message(STATUS "Found DLL ${libdir}/${librootname}.dll")
            endif ()
            set(${scipkgreg}_DLLS ${${scipkgreg}_DLLS} ${libdir}/${librootname}.dll)
          else ()
            if (DEBUG_CMAKE)
              message(STATUS "${librootname} has no accompanying dll.")
            endif ()
          endif ()
        endif ()
      endif ()
    endforeach ()
  endif ()

####################################################################
#
# Determine if this package should be marked as FOUND or not
# At the moment, we mark it as found if ANYTHING was found
# Per http://www.cmake.org/Wiki/scimake:How_To_Find_Installed_Software,
# The convention is to capitalize the _FOUND variable.
#
####################################################################

  if (${scipkgreg}_FOUND_SOME_LIBRARY OR ${scipkgreg}_FOUND_SOME_HEADER OR
       ${scipkgreg}_FOUND_SOME_EXECUTABLE OR ${scipkgreg}_FOUND_SOME_DLL OR
       ${scipkgreg}_FOUND_SOME_FILE)
    set(${scipkguc}_FOUND TRUE)
    if (DEBUG_CMAKE OR NOT ${scipkgreg}_FIND_QUIETLY)
      message(STATUS "Found ${scipkgreg}.")
      SciPrintCMakeResults(${scipkgreg})
    endif ()
  else ()
    set(${scipkguc}_FOUND FALSE)
    if (DEBUG_CMAKE)
      message(STATUS "Failed to find package ${scipkgreg}")
      SciPrintCMakeResults(${scipkgreg})
    endif ()
# If this was marked as a required package, fail
    if (${scipkgreg}_FIND_REQUIRED)
      message(FATAL_ERROR "Unable to find required package ${scipkgreg} - failing")
    endif ()
  endif ()

  message("--------- SciFindPackage done with ${TFP_PACKAGE} -----------")
  message("")
  # if (DEBUG_CMAKE)
    message(STATUS "${scipkguc}_FOUND = ${${scipkguc}_FOUND}.")
  # endif ()

endmacro(SciFindPackage)

