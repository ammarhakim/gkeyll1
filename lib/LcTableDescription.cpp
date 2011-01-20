/**
 * @file	LcTableDescription.cpp
 *
 * @brief	Description of a single value in a Lua table.
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
#include <LcTableDescription.h>

namespace Lucee
{
  TableDescription::TableDescription(const std::string& nm)
    : name(nm)
  {
  }

  void
  TableDescription::checkAndSet(Lucee::LuaTable& tbl)
  {
// check and set values
    checkAndSetValues<int>(tbl);
    checkAndSetValues<double>(tbl);
    checkAndSetValues<std::string>(tbl);
// check and set vectors
    checkAndSetVectors<int>(tbl);
    checkAndSetVectors<double>(tbl);
    checkAndSetVectors<std::string>(tbl);
  }
}
