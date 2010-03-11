/**
 * @file	LcRteRegistry.h
 *
 * @brief	Class for registering RTE solver object.
 *
 * @version	$Id: LcRteRegistry.h 321 2010-03-04 23:13:16Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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

