/**
 * @file	LcGenericFactory.h
 *
 * @brief	Templated class to serve as base class as a object 
 *                 factory.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_GENERIC_FACTORY_H
#define LC_GENERIC_FACTORY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee include
#include <LcLuaTable.h>
#include <LcSolverIfc.h>

namespace Lucee
{
/**
 * This class serves as a "proxy factory" to create objects based on
 * input files for those classes that can not initialize themselves
 * from the readInput() method. This is generally the case with most
 * classes as they have default constructors that need to be called
 * instead or before calling the readInput() method.
 */
  template <class BASEOBJ>
  class GenericFactory
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/** Destructor */
      virtual ~GenericFactory()
      {}

/**
 * Bootstrap method: Read input from specified table.
 *
 * @param tbl Table of input values.
 */
      virtual void readInput(Lucee::LuaTable& tbl) = 0;

/**
 * Create a new object derived from BASEOBJ and return pointer to
 * created object.
 *
 * @param solver Reference to containing solver.
 * @return pointer to class derived from BASEOBJ.
 */
      virtual BASEOBJ* create(const Lucee::SolverIfc& solver) = 0;
  };
}

#endif // LC_GENERIC_FACTORY_H
