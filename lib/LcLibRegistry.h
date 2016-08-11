/**
 * @file	LcLibRegistry.h
 *
 * @brief	Method for registering basic library object.
 */

#ifndef LC_LIB_REGISTRY_H
#define LC_LIB_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register basic solver objects and modules.
 */
  void registerLibObjects(Lucee::LuaState& L);
}

#endif // LC_LIB_REGISTRY_H
