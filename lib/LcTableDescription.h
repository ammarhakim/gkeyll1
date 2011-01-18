/**
 * @file	LcTableDescription.h
 *
 * @brief	Description of a single value in a Lua table.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_TABLE_DESCRIPTION_H
#define LC_TABLE_DESCRIPTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcValueDescription.h>

// std includes
#include <map>
#include <string>

namespace Lucee
{
/**
 * This class describes a Lua table that can appear in a Lua
 * script. It mainly served as a validation mechanism for the Lua
 * table and also for autmatically generating documentation for each
 * table.
 */
  class TableDescription
  {
    public:
/**
 * Create a new table description object.
 *
 * @param nm Name of table.
 */
      TableDescription(const std::string& nm);

    private:
/** Name of table */
      std::string name;
  };
}

#endif // LC_TABLE_DESCRIPTION_H
