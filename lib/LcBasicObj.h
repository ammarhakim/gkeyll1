/**
 * @file	LcBasicObj.h
 *
 * @brief	Interface class for basic Lucee objects.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_BASIC_OBJ_H
#define LC_BASIC_OBJ_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcLuaFuncMap.h>
#include <LcLuaTable.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Represents a simple object in Lucee that can be initialized from a
 * Lua script. Its main purpose is to provide a method of holding
 * pointers to derived type objects in maps. It also provides a
 * default implementation of a static function to add Lua callable
 * methods to Lucee.
 */
  class BasicObj
  {
    public:
/**
 * Create an object.
 */
      BasicObj();

/**
 * Create an object with specified name.
 *
 * @param nm Name of object.
 */
      BasicObj(const std::string& nm);

/**
 * Destroy object.
 */
      virtual ~BasicObj();

/**
 * Default method that performs registration of Lua functions. This
 * function does nothing: if derived classes need to register Lua
 * callable functions they must provide this method. Methods should be
 * added to the lfm object by calling the appendFunc() method.
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm);

/**
 * Get name of solver.
 *
 * @return Name of solver.
 */
      std::string getName() const;

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl) = 0;

    protected:
/**
 * Set name of object. This should be called by derived classes to set
 * their names.
 *
 * @param nm Name of solver.
 */
      void setName(const std::string& nm);
      
    private:
/** Name of the object */      
      std::string nm;
  };
}

#endif // LC_BASIC_OBJ_H
