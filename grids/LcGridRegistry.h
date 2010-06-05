/**
 * @file	LcGridRegistry.h
 *
 * @brief	Class for registering grid objects.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_GRID_REGISTRY_H
#define LC_GRID_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register grid specific objects and modules.
 */
  void registerGridObjects(Lucee::LuaState& L);
}

#endif // LC_GRID_REGISTRY_H
