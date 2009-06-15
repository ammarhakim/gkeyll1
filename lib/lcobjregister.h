/**
 * @file	lcobjregister.h
 *
 * @brief	Class for handling object registration
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_OBJ_REGISTER_H
#define LC_OBJ_REGISTER_H

// lib includes
#include <lcexcept.h>

// etc includes
#include <loki/Factory.h>
#include <loki/Singleton.h>

namespace Lucee
{
/**
 * Class to perform registration.
 */
  template<class D, class B>
  class ObjRegister
  {
    public:
/**
 * Register a new object to be created by its given name.
 *
 * @param nm Name by which object will be created.
 */
      ObjRegister(const std::string& nm)
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

#endif // LC_OBJ_REGISTER_H
