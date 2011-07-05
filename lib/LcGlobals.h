/**
 * @file	LcGlobals.h
 *
 * @brief	Class to hold global data used in various parts of Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2011, Ammar Hakim.
 */

#ifndef LC_GLOBALS_H
#define LC_GLOBALS_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to hold global data used in Lucee
 */
  struct Globals
  {
/** Output prefix for files */
      std::string outPrefix;
  };
}

#endif // LC_GLOBALS_H
