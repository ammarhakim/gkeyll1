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
 * Class to represent a Lua table. This constructor assumes that the
 * Lua table is on top of the stack.
 */
  class LuaTable
  {
    public:
/**
 * Create a new table object.
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
 * Get the type of this table. Not all tables will have a type field.
 *
 * @return type of this table.
 */
      std::string getType();

/**
 * Get the kind of this table. Not all tables will have a kind field.
 *
 * @return kind of this table.
 */
      std::string getKind();

/**
 * Get all numbers in table.
 *
 * @return numbers in table.
 */
      std::vector<double> getAllNumbers();

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

/**
 * Get a table inside this table.
 * 
 * @param nm Name of table to fetch.
 * @return table object.
 */
      LuaTable getTable(const std::string& nm);

/**
 * Check if string is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise
 */
      bool hasString(const std::string& key);

/**
 * Check if number is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasNumber(const std::string& key);

/**
 * Check if vector of strings is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasStrVec(const std::string& key);

/**
 * Check if vector of numbers is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasNumVec(const std::string& key);

/**
 * Check if a table is inside this table.
 * 
 * @param nm Name of table to fetch.
 * @return true if exists, false otherwise.
 */
      bool hasTable(const std::string& nm);

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
