/**
 * @file	LcLuaFunction.h
 *
 * @brief	Class to represent Lua function.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_FUNCTION_H
#define LC_LUA_FUNCTION_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaState.h>

// lua includes
#include <lua.hpp>

// std includes
#include <string>

namespace Lucee
{
/**
 * Class to represent a Lua function. This constructor assumes that the
 * Lua function is on top of the stack.
 */
  class LuaFunction
  {
    public:
/**
 * Create a new function object.
 *
 * @param L Lua state.
 * @param nm Name of function.
 */
      LuaFunction(Lucee::LuaState& L, const std::string& nm);

/**
 * Destroy table: frees Lua table and allows the garbage collector to
 * free memory.
 */
      ~LuaFunction();

/**
 * Evaluate function and return result.
 *
 * @param t Time to evaluate function.
 * @param xyz Coordinates of location.
 * @return value at (t, xyz)
 */
      double eval(double t, double xyz[3]);

    private:
/** Reference to lua state */
      Lucee::LuaState& L;
/** Name of the function */
      std::string name;
/** Pointer to Lua table */
      int ref;
  };
}

#endif // LC_LUA_TABLE_H
