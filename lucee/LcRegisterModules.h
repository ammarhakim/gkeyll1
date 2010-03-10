/**
 * @file	LcRegisterModules.h
 *
 * @brief	Class for registering all modules.
 *
 * @version	$Id: LcRegisterModules.h 321 2010-03-04 23:13:16Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_REGISTER_MODULES_H
#define LC_REGISTER_MODULES_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcRteRegistry.h>

namespace Lucee
{
  void registerModules(Lucee::LuaState& L);
}

#endif // LC_REGISTER_MODULES_H
