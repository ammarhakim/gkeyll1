/**
 * @file	LcLuceeObj.cpp
 *
 * @brief	Base class for generic Lua-callable objects in Lucee.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuceeObj.h>

namespace Lucee
{
// set module name
  const char *LuceeObj::id = "Lucee";
}
