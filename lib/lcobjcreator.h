/**
 * @file	lcobjcreator.h
 *
 * @brief	Class for handling object creation.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_OBJ_CREATOR_H
#define LC_OBJ_CREATOR_H

// lib includes
#include <lcexcept.h>

// etc includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

namespace Lucee
{
/**
 * Class to create objects. This class is used in conjunction with the
 * Lucee::ObjRegister class. Once an object has been registered this
 * class can be used to create derived classes using their names.
 */
  template <class B>
  class ObjCreator
  {
    public:
/**
 * Get a new object whose creator has the given name. The returned
 * object points to the base class.
 *
 * @param nm Name of the creator.
 * @return pointer to newly created object.
 */
      static B* getNew(const std::string& nm) 
      {
        return Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().CreateObject(nm);
      }

/**
 * Get a list of registered names.
 *
 * @return List of registered names.
 */
      static std::vector<std::string> registeredNames() 
      {
        return Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().RegisteredIds();
      }

/**
 * Unregister creator with given name.
 *
 * @param nm Name of creator to unregister.
 * @return bool true if unregister worked, false otherwise.
 */
      static bool unregister(const std::string& nm)
      {
        return Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().Unregister(nm);
      }
  };
}

#endif // LC_OBJ_CREATOR_H
