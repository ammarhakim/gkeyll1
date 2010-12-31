# $Id: TxSvnInfo.cmake 146 2010-10-21 09:37:02Z cary $
#
# For getting the svn revision of a directory unix

MACRO(Subversion_GET_VERSION dir var1 var2)
  IF (EXISTS ${dir}/.svn)
    MESSAGE(STATUS "Executing ${SVNVERSION_BIN} ${dir}")
    EXECUTE_PROCESS(COMMAND ${SVNVERSION_BIN} ${dir}
      OUTPUT_VARIABLE ${var1}
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    EXECUTE_PROCESS(COMMAND ${SVN_BIN} info ${dir}
      OUTPUT_FILE ${dir}/svninfo.txt
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  ELSE (EXISTS ${dir}/.svn)
    IF (EXISTS ${dir}/svninfo.txt)
      FILE(READ ${dir}/svninfo.txt ${var1})
      STRING(REGEX REPLACE "^(.*\n)?Revision: ([^\n]+).*"
        "\\2" ${var1} "${${var1}}")
    ELSE (EXISTS ${dir}/svninfo.txt)
      SET(${var1} "unknown")
    ENDIF (EXISTS ${dir}/svninfo.txt)
  ENDIF (EXISTS ${dir}/.svn)
  IF (EXISTS ${dir}/svninfo.txt)
    FILE(READ ${dir}/svninfo.txt SVNINFO)
    STRING(REGEX REPLACE "^(.*\n)?URL: ([^\n]+).*"
      "\\2" ${var2} "${SVNINFO}")
  ELSE (EXISTS ${dir}/svninfo.txt)
    SET(${var2} "unknown")
  ENDIF (EXISTS ${dir}/svninfo.txt)
ENDMACRO(Subversion_GET_VERSION)

