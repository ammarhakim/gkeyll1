######################################################################
#
# FindSciDocutils
#
# $Id: FindSciDocutils.cmake 259 2013-04-10 19:10:45Z jdelamere $
#
# Copyright 2010-2013 Tech-X Corporation.
# Arbitrary redistribution allowed provided this copyright remains.
#
# See LICENSE file (EclipseLicense.txt) for conditions of use.
#
######################################################################
include(FindPackageHandleStandardArgs)

find_program(RST2HTML_EXECUTABLE NAMES rst2html rst2html.py)
find_package_handle_standard_args(Docutils DEFAULT_MSG RST2HTML_EXECUTABLE)

find_program(RST2LATEX_EXECUTABLE NAMES rst2latex rst2latex.py)
find_package_handle_standard_args(Docutils DEFAULT_MSG RST2LATEX_EXECUTABLE)
