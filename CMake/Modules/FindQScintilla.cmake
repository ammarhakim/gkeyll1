#################################################
#
# Find module for QScintilla includes and lib
#
# $Id: FindQScintilla.cmake 160 2010-11-23 19:46:17Z mkoch $
#
#################################################

FIND_PATH(QScintilla_INCLUDE_DIR 
  qsciscintilla.h
  PATHS ${QT_INCLUDE_DIR}
  PATH_SUFFIXES Qsci
)

FIND_LIBRARY(QScintilla_LIBRARY
  qscintilla2
  PATHS ${QT_LIBRARY_DIR}
)

IF(QScintilla_INCLUDE_DIR AND QScintilla_LIBRARY)
  MESSAGE(STATUS "Found QScintilla")
  MESSAGE(STATUS "QScintilla_INCLUDE_DIR = ${QScintilla_INCLUDE_DIR}")
  MESSAGE(STATUS "QScintilla_LIBRARY = ${QScintilla_LIBRARY}")
ELSE(QScintilla_INCLUDE_DIR AND QScintilla_LIBRARY)
  IF(QScintilla_FIND_REQUIRED)
    MESSAGE(FATAL "Could not find QScintilla, make sure it has been installed to Qt directories")
  ENDIF(QScintilla_FIND_REQUIRED)
ENDIF(QScintilla_INCLUDE_DIR AND QScintilla_LIBRARY)
