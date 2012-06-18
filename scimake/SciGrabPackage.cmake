######################################################################
#
# SciInstallLibs:
#    Macros for automating the grabbing of packages installed by bilder
#    at the CMAKE_INSTALL_DIR/.. directory
#
# $Id: SciGrabPackage.cmake 1344 2012-05-11 18:10:13Z kruger $
#
# Copyright 2010-2012 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
######################################################################

include(CMakeParseArguments)

# SciGrabPackage.cmake
#  Consider hdf5.  SciFindPackage will normally look for either serial
#   or parallel, but want to grab: ser,par,cc4py at one go.  This
#   enables one to specify one of them, and then grab all of them.
#
# Args:
#   PACKAGE:  Basename of package.  (e.g., hdf5)
#   VERSIONS: List of versions to grab (e.g., ser,par,cc4py) 
#     Note that ser maps to search for PACKAGE, but all others
#     map to a search for PACKAGE-VERSION  (e.g., search for hdf5-par)
#   PKGDIR: Where one of the packages versions are installed 
#       (e.g., #   /contrib/hdf5-par)
#   INSTALLDIR: Where to install the packages
#   COMMON_VERSION: Of all the versions, which one goes into common area
#   COMMON_INSTALL: If invoked, put static libraries and binaries in
#        common version as well.
#
macro(SciGrabPackage)

# Parse out the args
  set(opts DEBUG;COMMON_INSTALL) # no-value args
  set(oneValArgs PACKAGE;PKGDIR;INSTALLDIR;COMMON_VERSION)
  set(multValArgs VERSIONS) # e.g., lists
  cmake_parse_arguments(FD "${opts}" "${oneValArgs}" "${multValArgs}" ${ARGN})

  ###
  ##  Basic idiot checks
  #
  if(NOT DEFINED FD_PACKAGE)
     message(WARNING "SciGrab called without defining package")
     return()
  endif()
  if(NOT DEFINED FD_PKGDIR)
     message(WARNING "SciGrab called without defining where to look for package")
     return()
  endif()
  if(NOT DEFINED FD_VERSIONS)
     message(WARNING "SciGrab called without defining versions to look for ")
     return()
  endif()
  if(NOT DEFINED FD_INSTALLDIR)
     message(WARNING "SciGrab called without defining installdir")
     return()
  endif()
  if (FD_DEBUG)
    message(STATUS "[SciGrabPackage]: PACKAGE= ${FD_PACKAGE} ")
    message(STATUS "[SciGrabPackage]: PKGDIR= ${FD_PKGDIR} ")
    message(STATUS "[SciGrabPackage]: VERSIONS= ${FD_VERSIONS} ")
    message(STATUS "[SciGrabPackage]: INSTALLDIR= ${FD_INSTALLDIR} ")
  endif()

  message("")
  message("--------- SciGrabPackage grabbing ${FD_PACKAGE} ---------")
  ###
  ##  Grab all installed versions (ser,par,...) of the installations
  ##  First figure out the basedir and see if we have all of the versions
  ##  To figure out the "common_version" is complicated because we have
  ##   to have the resolved name (hdf5 -> hdf5-<version>-ser)
  #
  get_filename_component(bilderdir ${FD_PKGDIR}/.. REALPATH)
  if(NOT DEFINED FD_COMMON_VERSION)
    set(FD_COMMON_VERSION "ser")
  endif()

  set (foundversions)
  foreach (version ${FD_VERSIONS})
     if (${version} STREQUAL "ser")
        set(packagedir ${bilderdir}/${FD_PACKAGE})
     else()
        set(packagedir ${bilderdir}/${FD_PACKAGE}-${version})
     endif()
     if (NOT EXISTS ${packagedir})
        message(STATUS "SciGrab cannot find ${FD_PACKAGE} at ${packagedir}")
     else()
        # Change to a full path
        get_filename_component(fullpkgdir ${packagedir} REALPATH)
        if (${FD_COMMON_VERSION} STREQUAL "${version}")
          get_filename_component(COMMON_VERSION_NAME ${fullpkgdir} NAME)
        endif()
        list(APPEND foundversions ${fullpkgdir})
        message(STATUS ${fullpkgdir})
     endif()
  endforeach()
  if (NOT foundversions )
     message(STATUS "SciGrab could not find ${PACKAGE}-(${VERSIONS}) in ${bilderdir}")
     return()
  endif()

  ###
  ## Once found, now we can install it
  ## Here are the rules:
  ##  if verdir == COMMON_VERSION
  ##     Shared libs: Always install in ${FD_INSTALLDIR}/lib
  ##     if COMMON_INSTALL
  ##         Static libs and binaries: ${FD_INSTALLDIR}/{lib,bin} 
  ##     else
  ##         Static libs and binaries: ${FD_INSTALLDIR}/verdirname/{lib,bin}
  ##  else
  ##     Everything: ${FD_INSTALLDIR}/verdirname/{lib,bin}
  ## 
  #
  foreach (verdir ${foundversions})
    message(STATUS "DIR TO INSTALL" ${verdir})
    ###
    ##  Everything in its right place
    ##  as described above
    #
    get_filename_component(verdirname ${verdir} NAME)
    message(STATUS "TMP COMMON_VERSION_NAME:: ${COMMON_VERSION_NAME}")
    message(STATUS "TMP:: verdername  ${verdirname}")
    message(STATUS "TMP:: verdir  ${verdir}")
    if ("${verdirname}" STREQUAL "${COMMON_VERSION_NAME}")
      set(INSTSHLIBDIR   ${FD_INSTALLDIR}/lib)
      if(DEFINED FD_COMMON_INSTALL)
        set(INSTSTLIBDIR ${FD_INSTALLDIR}/lib)
        set(INSTBINDIR   ${FD_INSTALLDIR}/bin)
        set(INSTINCDIR   ${FD_INSTALLDIR})
      endif()
      # This is for non bin or lib directories
      set(INSTDIRGEN   ${FD_INSTALLDIR}/${verdirname})
      #set(INSTDIRGEN   ${FD_INSTALLDIR})
      set(specialdirs "bin;lib")
    else()
      set(INSTSHLIBDIR ${FD_INSTALLDIR}/${verdirname}/lib)
      set(INSTSTLIBDIR ${FD_INSTALLDIR}/${verdirname}/lib)
      set(INSTBINDIR   ${FD_INSTALLDIR}/${verdirname}/bin)
      set(INSTINCDIR   ${FD_INSTALLDIR}/${verdirname})
      set(INSTDIRGEN   ${FD_INSTALLDIR}/${verdirname})
      set(specialdirs "bin;lib")
    endif()
    ###
    ## Install the bins and the libraries
    #
    install(CODE
      "include(${SCICMAKE_DIR}/SciFixDyLibEx.cmake)
      SciFixDyLibEx(FROM_DIR ${verdir} FROM_BINDIR bin FROM_LIBDIR lib 
           BINDIR ${INSTBINDIR} LIBDIR ${INSTSTLIBDIR} SHLIBDIR ${INSTSHLIBDIR} 
           INSTALL_PYTHON ${ENABLE_PYTHON} DEBUG)"
    )
    ###
    ## Install everything else
    #
    file(GLOB vsubdirs ${verdir}/*)
    if (FD_DEBUG)
      message(STATUS "[SciGrabPackage]: vsubdirs= ${vsubdirs} ")
    endif ()
    foreach (vsdir ${vsubdirs})
      install(CODE "MESSAGE(STATUS \"SciGrabPackage: ${vsdir}\")")
      message(STATUS "[SciGrabPackage]: vsdir= ${vsdir} ")
      get_filename_component(vsdirname ${vsdir} NAME_WE)
      if ("${vsdirname}" STREQUAL "bin" OR "${vsdirname}" STREQUAL "lib")
        if (FD_DEBUG)
          message(STATUS "[SciGrabPackage]: found bin or lib= ${vsdirname} ")
        endif ()
      else()
        if ("${vsdirname}" STREQUAL "include")
          install(DIRECTORY ${vsdir} DESTINATION ${INSTINCDIR} USE_SOURCE_PERMISSIONS)
        else ()
          if(IS_DIRECTORY ${vsdir})
            install(DIRECTORY ${vsdir} DESTINATION ${INSTDIRGEN} USE_SOURCE_PERMISSIONS)
          else()
            message(STATUS "[SciGrabPackage]: Not installing ${vsdir} ")
          endif()
        endif()
      endif ()
    endforeach()

  endforeach()
endmacro()

