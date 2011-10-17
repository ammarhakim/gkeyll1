/**
 * @file	LcLuceeMod.cpp
 *
 * @brief	Base class for generic Lua-callable modects in Lucee.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuceeMod.h>

namespace Lucee
{
// set module name
  const char *LuceeMod::id = "Lucee";
}
