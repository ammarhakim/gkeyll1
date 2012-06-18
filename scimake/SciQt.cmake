######################################################################
#
# Shared CMakeLists.txt for all Composers
#
# $Id: SciQt.cmake 1324 2012-04-21 13:31:24Z cary $
#
#####################################################################

# It would be nice to find the Qt libraries directly in VisIt, but
# FindQt4 relies on qmake, which VisIt does not install.
# So we find them here, and then we replace with the VisIt versions
# if they are found.
set(QT_REQUIRED_LIBRARIES
  QtCore QtGui QtOpenGL QtXml QtXmlPatterns QtNetwork QtWebKit
)

# The following is included here for the call to QT4_ADD_RESOURCES
find_package(SciQt4 COMPONENTS ${QT_REQUIRED_LIBRARIES})

# Fix up the libraries on Linux to remove debug versions and limiters.
# We use these libraries even if using VisIt, as later we will change
# these to the libraries installed with VisIt.
if (0) # Done by FindSciQt4?
if (NOT APPLE)
  set(tmp ${QT_LIBRARIES})
  set(QT_LIBRARIES)
  foreach (qtlib ${tmp})
      if (NOT ${qtlib} STREQUAL optimized AND NOT ${qtlib} MATCHES debug)
        set(QT_LIBRARIES ${QT_LIBRARIES} ${qtlib})
      endif ()
  endforeach ()
endif ()
endif ()

# Add in phonon if found in same directory
if (QT_PHONON_LIBRARY_RELEASE)
  get_filename_component(phonondir ${QT_PHONON_LIBRARY_RELEASE}/.. REALPATH)
  get_filename_component(qtcoredir ${QT_QTCORE_LIBRARY_RELEASE}/.. REALPATH)
  if (${phonondir} STREQUAL ${qtcoredir})
    message(STATUS "[SciComposerBase]: Adding ${QT_PHONON_LIBRARY_RELEASE} to QT_LIBRARIES.")
    set(QT_LIBRARIES ${QT_LIBRARIES} ${QT_PHONON_LIBRARY_RELEASE})
    if (WIN32)
      get_filename_component(phonondll ${phonondir}/../bin/phonon4.dll REALPATH)
      if (EXISTS ${phonondll})
        message(STATUS "[SciComposerBase]: Adding ${phonondll} to QT_DLLS.")
        set(QT_DLLS ${QT_DLLS} ${phonondll})
      else ()
        message(STATUS "[SciComposerBase]: phonon4.dll not found in ${phonondir}/../bin.")
      endif ()
    endif ()
  else ()
    message(STATUS "[SciComposerBase]: Not adding ${QT_PHONON_LIBRARY_RELEASE} to QT_LIBRARIES since not from same installation.")
  endif ()
endif ()

# Here's a trick - we want to pass the list of qt libraries to
# our postbuild script.  But it's a semi-colon delimited list
# and we need space-delimited.  So convert
set(QT_REQUIRED_LIBRARIES_STRING CACHE STRING "List of required Qt libraries")
foreach (arg ${QT_REQUIRED_LIBRARIES})
  set(QT_REQUIRED_LIBRARIES_STRING "${QT_REQUIRED_LIBRARIES_STRING} ${arg}")
endforeach ()
# JRC: could just do a substitution

# If not using VisIt, add includes from the Qt installation.
if (NOT ENABLE_VISIT)

# Include here unless looking for VisIt, in which case we will change
  include_directories(
    ${QT_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${QT_QSCIML_INCLUDE_DIR}
    ${QT_QTNETWORK_INCLUDE_DIR}
    ${QT_QTWEBKIT_INCLUDE_DIR}
    ${QT_QTGUI_INCLUDE_DIR}
    ${QT_QTOPENGL_INCLUDE_DIR}
  )
  if (QT_PHONON_INCLUDE_DIR)
    include_directories(${QT_PHONON_INCLUDE_DIR})
  endif ()
endif ()

macro(SciTranspPreprocess)
# Install the Qt plugins
set(QT_PLUGIN_INSTDIR)
if (APPLE)
  set(QT_PLUGIN_INSTDIR "${FLAVOR_INSTALL_ROOT}")
else ()
  set(QT_PLUGIN_INSTDIR "${FLAVOR_INSTALL_BINDIR}")
endif ()
install(CODE "MESSAGE(STATUS \"******** Installing the Qt plugins into ${QT_PLUGIN_INSTDIR}\")")
install(DIRECTORY "${QT_PLUGINS_DIR}"
  DESTINATION "${QT_PLUGIN_INSTDIR}"
  PATTERN "*_debug.dylib" EXCLUDE
  PATTERN "*.lib" EXCLUDE
  PATTERN "*d4.*" EXCLUDE
  PATTERN "accessible" EXCLUDE
  PATTERN "bearer" EXCLUDE
  PATTERN "codecs" EXCLUDE
  PATTERN "qmltooling" EXCLUDE
  PATTERN "libqdeclarativeview.dylib" EXCLUDE
  PATTERN "libphononwidgets.dylib" EXCLUDE
  PATTERN "graphicssystems" EXCLUDE
  PATTERN "phonon_backend" EXCLUDE
  PATTERN "*svg.*" EXCLUDE
  PERMISSIONS OWNER_WRITE OWNER_READ
              GROUP_WRITE GROUP_READ
              WORLD_READ
)

# Fix executable and visit libs references on OS X
if (APPLE AND NOT ENABLE_VISIT)

  if (DEBUG_CMAKE)
    message(STATUS "[SciComposerBase]: QT_LIBRARY_DIR = ${QT_LIBRARY_DIR}.")
    install(CODE
      "
      include(${CMAKE_SOURCE_DIR}/scimake/SciFixDylibs.cmake)
      # SciFixDylibs(TOFIXOBJ ${APPLICATION_NAME} INSTALL_SUBDIR ${FLAVOR_INSTALL_BINDIR} NEW_LIBDIR @loader_path/../VisIt/${VisIt_ARCH_SUBDIR}/lib DEBUG)
      # SciFixDylibs(INSTALL_SUBDIR ${FLAVOR_INSTALL_ROOT}/plugins/imageformats OLD_LIBDIR ${QT_LIBRARY_DIR} NEW_LIBDIR @loader_path/../../VisIt/${VisIt_ARCH_SUBDIR}/lib DEBUG)
      # SciFixDylibs(INSTALL_SUBDIR ${FLAVOR_INSTALL_ROOT}/plugins/designer OLD_LIBDIR ${QT_LIBRARY_DIR} NEW_LIBDIR @loader_path/../../VisIt/${VisIt_ARCH_SUBDIR}/lib DEBUG)
      "
    )
  else ()
    install(CODE
      "
      include(${CMAKE_SOURCE_DIR}/scimake/SciFixDylibs.cmake)
      # SciFixDylibs(TOFIXOBJ ${APPLICATION_NAME} INSTALL_SUBDIR ${FLAVOR_INSTALL_BINDIR} NEW_LIBDIR @loader_path/../VisIt/${VisIt_ARCH_SUBDIR}/lib)
      # SciFixDylibs(INSTALL_SUBDIR ${FLAVOR_INSTALL_ROOT}/plugins/imageformats OLD_LIBDIR ${QT_LIBRARY_DIR} NEW_LIBDIR @loader_path/../../VisIt/${VisIt_ARCH_SUBDIR}/lib)
      # SciFixDylibs(INSTALL_SUBDIR ${FLAVOR_INSTALL_ROOT}/plugins/designer OLD_LIBDIR ${QT_LIBRARY_DIR} NEW_LIBDIR @loader_path/../../VisIt/${VisIt_ARCH_SUBDIR}/lib)
      "
    )
  endif ()

endif ()
# Qt.conf(tells the app where to load Qt plugins)
# Needed on ALL PLATFORMS
if (APPLE)
  install(CODE "MESSAGE(STATUS \"******** Installing qt.conf in ${FLAVOR_INSTALL_ROOT}/Resources\")")
  install(FILES ${COMPOSERTOOLKIT_DIR}/bin/qt.conf
          DESTINATION ${FLAVOR_INSTALL_ROOT}/Resources
          PERMISSIONS OWNER_WRITE OWNER_READ
                      GROUP_WRITE GROUP_READ
                      WORLD_READ
  )
else ()
  install(CODE "MESSAGE(STATUS \"******** Installing qt.conf in ${FLAVOR_INSTALL_BINDIR}\")")
  install(FILES ${COMPOSERTOOLKIT_DIR}/bin/qt.conf
          DESTINATION ${FLAVOR_INSTALL_BINDIR}
          PERMISSIONS OWNER_WRITE OWNER_READ
                      GROUP_WRITE GROUP_READ
                      WORLD_READ
  )
endif ()

endmacro()

