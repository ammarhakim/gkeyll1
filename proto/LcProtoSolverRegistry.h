/**
 * @file	LcProtoSolverRegistry.h
 *
 * @brief	Method for registering basic solver object.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObjRegistry.h>

namespace Lucee
{
  void
  registerProtoSolverObjects(Lucee::LuaState& L);
}
