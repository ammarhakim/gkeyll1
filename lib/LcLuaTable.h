/**
 * @file	LcLuaTable.h
 *
 * @brief	Class to represent Lua table.
 */

#ifndef LC_LUA_TABLE_H
#define LC_LUA_TABLE_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcLuaState.h>
#include <LcPointerHolder.h>

// lua includes
#include <lua.hpp>

// std include
#include <map>
#include <string>
#include <vector>

// #include <iostream>

// #define SHOW_LUA_STACK_SIZE(nm, L) \
//     do { \
//       std::cout << "Stack size in " << nm << " is " << lua_gettop(L) << std::endl; \
//     } while (0);

// #define SHOW_LUA_STACK_SIZE2(L) \
//     do { \
//     std::cout << "End stack size is " << lua_gettop(L) << std::endl; \
//     } while (0);

/** Debug macro for tracking Lua stack */
#define SHOW_LUA_STACK_SIZE(nm, L) \
    do {} while (0);
/** Debug macro for tracking Lua stack */
#define SHOW_LUA_STACK_SIZE2(L) \
    do {} while (0);

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
 * Get a reference to the LuaState object connected to this table.
 *
 * @param reference to the LuaState object.
 */
      Lucee::LuaState& getLuaState();

/**
 * Get the type of this table. Not all tables will have a type field.
 *
 * @return type of this table.
 */
      std::string getType() const;

/**
 * Get the kind of this table. Not all tables will have a kind field.
 *
 * @return kind of this table.
 */
      std::string getKind() const;

/**
 * Get all numbers in table.
 *
 * @return numbers in table.
 */
      std::vector<double> getAllNumbers() const;

/**
 * Get all strings in table.
 *
 * @return strings in table.
 */
      std::vector<std::string> getAllStrings() const;

/**
 * Get all function references in table.
 *
 * @return list of function references in table
 */
      std::vector<int> getAllFunctionRefs() const;

/**
 * Get all objects in table as a vector of pointers.
 *
 * @return objects in table.
 */
      template <typename T>
      std::vector<T*> getAllObjects() const
      {
        SHOW_LUA_STACK_SIZE("getAllObjects", L);
        std::vector<T*> res;
// push table object on stack
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
        int t = lua_gettop(L);
        lua_pushnil(L);
        while (lua_next(L, t) != 0) 
        {
          if (lua_type(L, -1) == LUA_TUSERDATA)
          {
            Lucee::PointerHolder<T> *ph =
              (Lucee::PointerHolder<T>*) lua_touserdata(L, -1);
            res.push_back(ph->pointer);
          }
          lua_pop(L, 1);
        }
        lua_pop(L, 1);
        SHOW_LUA_STACK_SIZE2(L);
        return res;
      }

/**
 * Get a string from table.
 *
 * @param key Key in table.
 * @return string corresponding to key.
 */
      std::string getString(const std::string& key) const;

/**
 * Get a number from table.
 *
 * @param key Key in table.
 * @return number corresponding to key.
 */
      double getNumber(const std::string& key) const;

/**
 * Get a boolean from table.
 *
 * @param key Key in table.
 * @return true or false
 */
      bool getBool(const std::string& key) const;

/**
 * Get user data (i.e. pointer to Lucee object) from table. The object
 * type must be specified as the template parameter.
 *
 * @param key Key in table.
 * @return reference to object.
 */
      template <typename T>
      T& getObject(const std::string& key) const
      {
// push table object on stack
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
        lua_pushstring(L, key.c_str());
// get data from table
        lua_gettable(L, -2);
        if (lua_type(L, -1) != LUA_TUSERDATA)
        {
          lua_pop(L, 2);
          Lucee::Except lce("LuaTable::getObject: ");
          lce << key << " is not a Lucee object";
          throw lce;
        }
        Lucee::PointerHolder<T> *ph =
          (Lucee::PointerHolder<T>*) lua_touserdata(L, -1);
        
        lua_pop(L, 2);
        SHOW_LUA_STACK_SIZE2(L);
        return *ph->pointer;
      }

/**
 * Get user data (i.e. pointer to Lucee object) from table. The object
 * type must be specified as the template parameter. This method
 * checks the type of the object, making sure it corresponds to the
 * base class for the class system.
 *
 * @param key Key in table.
 * @return reference to object.
 */
      template <typename T>
      T& getObjectAsBase(const std::string& key) const
      {
// push table object on stack
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
        lua_pushstring(L, key.c_str());
// get data from table
        lua_gettable(L, -2);
        if (lua_type(L, -1) != LUA_TUSERDATA)
        {
          lua_pop(L, 2);
          Lucee::Except lce("LuaTable::getObjectAsBase: ");
          lce << key << " is not a Lucee object";
          throw lce;
        }
        T *obj = Lucee::PointerHolder<T>::getObjAsBase(L, -1);
        
        lua_pop(L, 2);
        SHOW_LUA_STACK_SIZE2(L);
        return *obj;
      }

/**
 * Get user data (i.e. pointer to Lucee object) from table. The object
 * type must be specified as the template parameter. This method
 * checks the type of the object, making sure it corresponds to the
 * derived class for the class system.
 *
 * @param key Key in table.
 * @return reference to object.
 */
      template <typename T>
      T& getObjectAsDerived(const std::string& key) const
      {
// push table object on stack
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
        lua_pushstring(L, key.c_str());
// get data from table
        lua_gettable(L, -2);
        if (lua_type(L, -1) != LUA_TUSERDATA)
        {
          lua_pop(L, 2);
          Lucee::Except lce("LuaTable::getObjectAsDerived: ");
          lce << key << " is not a Lucee object";
          throw lce;
        }
        T *obj = Lucee::PointerHolder<T>::getObjAsDerived(L, -1);
        
        lua_pop(L, 2);
        SHOW_LUA_STACK_SIZE2(L);
        return *obj;
      }

/**
 * Get vector of strings from table.
 *
 * @param key Key in table.
 * @return vector of string corresponding to key.
 */
      std::vector<std::string> getStrVec(const std::string& key) const;

/**
 * Get vector of numbers from table.
 *
 * @param key Key in table.
 * @return vector of numbers corresponding to key.
 */
      std::vector<double> getNumVec(const std::string& key) const;

/**
 * Get vector of booleans from table.
 *
 * @param key Key in table.
 * @return vector of booleans corresponding to key.
 */
      std::vector<bool> getBoolVec(const std::string& key) const;

/**
 * Get a table inside this table.
 * 
 * @param nm Name of table to fetch.
 * @return table object.
 */
      LuaTable getTable(const std::string& nm) const;

/**
 * Get reference to specified Lua function.
 *
 * @param nm Name of function.
 * @return reference reference to function.
 */
      int getFunctionRef(const std::string& nm) const;

/**
 * Check if string is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise
 */
      bool hasString(const std::string& key) const;

/**
 * Check if number is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasNumber(const std::string& key) const;

/**
 * Check if boolean is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasBool(const std::string& key) const;

/**
 * Check if vector of strings is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasStrVec(const std::string& key) const;

/**
 * Check if vector of numbers is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasNumVec(const std::string& key) const;

/**
 * Check if vector of booleans is in table.
 *
 * @param key Key in table.
 * @return true if exists, false otherwise.
 */
      bool hasBoolVec(const std::string& key) const;

/**
 * Check if a table is inside this table.
 * 
 * @param nm Name of table to fetch.
 * @return true if exists, false otherwise.
 */
      bool hasTable(const std::string& nm) const;

/**
 * Check if a function is inside this table.
 * 
 * @param nm Name of function to check.
 * @return true if exists, false otherwise.
 */
      bool hasFunction(const std::string& nm) const;

/**
 * Check if a Lucee object exists in this table. The type of the
 * object must be passed as a template parameter.
 * 
 * @param key Name of object to check.
 * @return true if exisit, false otherwise.
 */
      template <typename T>
      bool hasObject(const std::string& key) const
      {
// push table object on stack
        lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
        lua_pushstring(L, key.c_str());
// get data from table
        lua_gettable(L, -2);
        if (lua_type(L, -1) != LUA_TUSERDATA)
        {
          lua_pop(L, 2);
          return false;
        }
        lua_pop(L, 2);
        SHOW_LUA_STACK_SIZE2(L);
        return true;
      }

/**
 * Get list of all table names with specified type.
 *
 * @param type Type of table.
 * @return list of table names.
 */
      std::vector<std::string> getNamesOfType(const std::string& type) const;

    private:
/** Reference to lua state */
      Lucee::LuaState& L;
/** Name of the table */
      std::string name;
/** Pointer to Lua table */
      int ref;
/** Map of "types" to table names */
      std::map<std::string, std::vector<std::string> > typeMap;

/**
 * Loop over all tables and create the typeMap.
 */
      void createTypeMap();

/**
 * Add a new variable to type map.
 */
      void addToTypeMap(const std::string& var, const std::string& type);
  };
}

#endif // LC_LUA_TABLE_H
