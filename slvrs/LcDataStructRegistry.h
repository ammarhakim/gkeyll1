/**
 * @file	LcDataStructRegistry.h
 *
 * @brief	Class for registering data-struct objects.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_DATA_STRUCT_REGISTRY_H
#define LC_DATA_STRUCT_REGISTRY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
/**
 * Register data-struct specific objects and modules.
 */
  void registerDataStructObjects(Lucee::LuaState& L);
}

#endif // LC_DATA_STRUCT_REGISTRY_H
