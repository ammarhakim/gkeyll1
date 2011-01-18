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
}
