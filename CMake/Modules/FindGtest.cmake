# $Id: FindGtest.cmake 165 2010-11-24 21:16:00Z amyx $

TxFindPackage("gtest" "" "gtest.h"
  "gtest" ""
  "include/gtest" "lib/${CXX_COMP_LIB_SUBDIR};lib")

