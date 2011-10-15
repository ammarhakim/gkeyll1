/**
 * @file	LcLuceeObj.h
 *
 * @brief	Base class for generic Lua-callable objects in Lucee.
 */

#ifndef LC_LUCEE_OBJ_H
#define LC_LUCEE_OBJ_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcBasicObj.h>

namespace Lucee
{
/**
 * Base class for generic Lua-callable objects in Lucee. This class
 * essentially provides the module "Lucee" and allows namespacing Lua
 * creatable objects.
 */
  class LuceeObj : Lucee::BasicObj
  {
    public:
/** Class id: this is used by registration system */
      static const char *id;
  };
}

#endif // LC_LUCEE_OBJ_H
