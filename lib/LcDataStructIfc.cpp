/**
 * @file	LcDataStructIfc.cpp
 *
 * @brief	Base class for all data in Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcDataStructIfc.h>

namespace Lucee
{
// set madule name
  const char *DataStructIfc::id = "DataStruct";

  DataStructIfc::DataStructIfc()
    : Lucee::BasicObj("DataStruct")
  {
  }

  DataStructIfc::~DataStructIfc()
  {
  }

  void
  DataStructIfc::readInput(Lucee::LuaTable& tbl)
  {
  }

  DataStructIfc*
  DataStructIfc::clone() const
  {
    return 0;
  }
}
