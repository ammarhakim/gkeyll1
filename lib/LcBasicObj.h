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
#include <LcLuaTable.h>

// std includes
#include <string>

namespace Lucee
{
/**
 * Represents the simplest possible object in Lucee. This object can
 * be initialized from a lua table. Its main purpose is to provide a
 * method of holding pointers to derived type objects in maps.
 */
  class BasicObj
  {
    public:
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
