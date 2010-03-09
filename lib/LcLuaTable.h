/**
 * @file	LcLuaTable.h
 *
 * @brief	Class to represent Lua table.
 *
 * @version	$Id: LcLuaTable.h 157 2009-08-26 17:27:55Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_LUA_TABLE_H
#define LC_LUA_TABLE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaState.h>

// lua includes
#include <lua.hpp>

// std include
#include <string>
#include <vector>

namespace Lucee
{
/**
 * Class to represent a Lua table.
 */
  class LuaTable
  {
    public:
/**
 * Create a new table object. This constructor assumes that the Lua
 * table is on top of the stack.
 *
 * @param L Lua state.
 * @param nm Name of table.
 */
      LuaTable(Lucee::LuaState& L, const std::string& nm);

/**
 * Destroy table: frees Lua table and allows the garbage collector to
 * free memory.
 */
      ~LuaTable();

/**
 * Get a string from table.
 *
 * @param key Key in table.
 * @return string corresponding to key.
 */
      std::string getString(const std::string& key);

/**
 * Get a number from table.
 *
 * @param key Key in table.
 * @return number corresponding to key.
 */
      double getNumber(const std::string& key);

/**
 * Get vector of strings from table.
 *
 * @param key Key in table.
 * @return vector of string corresponding to key.
 */
      std::vector<std::string> getStrVec(const std::string& key);

/**
 * Get vector of numbers from table.
 *
 * @param key Key in table.
 * @return vector of numbers corresponding to key.
 */
      std::vector<double> getNumVec(const std::string& key);

    private:
/** Reference to lua state */
      Lucee::LuaState& L;
/** Name of the table */
      std::string name;
/** Pointer to Lua table */
      int ref;
  };
}

#endif // LC_LUA_TABLE_H
