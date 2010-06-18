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
#include <vector>

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
 * @param numOut size of result vector.
 */
      LuaFunction(Lucee::LuaState& L, const std::string& nm, unsigned numOut);

/**
 * Destroy table: frees Lua table and allows the garbage collector to
 * free memory.
 */
      ~LuaFunction();

/**
 * Evaluate function and return result.
 *
 * @param inp Vector of input values.
 * @return function output.
 */
      std::vector<double> eval(const std::vector<double>& inp);

    private:
/** Reference to lua state */
      Lucee::LuaState& L;
/** Name of the function */
      std::string name;
/** Pointer to Lua table */
      int ref;
/** Number of output results */
      unsigned numOut;
  };
}

#endif // LC_LUA_TABLE_H
