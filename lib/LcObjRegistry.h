/**
 * @file	LcObjRegistry.h
 *
 * @brief	Class for handling object registration
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_OBJ_REGISTRY_H
#define LC_OBJ_REGISTRY_H

// lucee includes
#include <LcExcept.h>

// loki includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

namespace Lucee
{
/**
 * Class to perform object registration.
 */
  template<class B, class D>
  class ObjRegistry
  {
    public:
/**
 * Register a new object to be created by its given name.
 *
 * @param nm Name by which object will be created.
 */
      ObjRegistry(const std::string& nm)
      {
        Loki::SingletonHolder<Loki::Factory<B, std::string> >
          ::Instance().Register(nm, getNew);
      }

    private:
/**
 * Return a newly allocated object. This is used as the creation
 * mechanism for the object.
 *
 * @return Newly allocated object.
 */
      static B* getNew() 
      {
        return new D;
      }
  };
}

#endif // LC_OBJ_REGISTRY_H
