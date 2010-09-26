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
#include <LcBasicObj.h>
#include <LcLuaTable.h>

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
  class GenericFactory : public Lucee::BasicObj
  {
    public:
/** Class id: this is used by the registration system */
      static const char *id;

/** Default construct */
      GenericFactory()
        : Lucee::BasicObj("name") 
      {
      }

/** Destructor */
      virtual ~GenericFactory()
      {
      }

/**
 * Method that performs registration of Lua functions. This function
 * simply redirects the call to the BASEOBJ class. 
 *
 * @param lfm Lua function map object.
 */
      static void appendLuaCallableMethods(Lucee::LuaFuncMap& lfm)
      {
        BASEOBJ::appendLuaCallableMethods(lfm);
      }

/**
 * Create a new object derived from BASEOBJ and return pointer to
 * created object.
 *
 * @return pointer to class derived from BASEOBJ.
 */
      virtual BASEOBJ* create() = 0;
  };
}

#endif // LC_GENERIC_FACTORY_H
