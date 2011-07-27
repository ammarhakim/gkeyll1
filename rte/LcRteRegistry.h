/**
 * @file	LcRteRegistry.h
 *
 * @brief	Class for registering RTE solver object.
 */

#ifndef LC_RTE_REGISTRY_H
#define LC_RTE_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register RTE specific objects and modules.
 */
  void registerRteObjects(Lucee::LuaState& L);
}

#endif // LC_RTE_REGISTRY_H

