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
#include <lcobjcreator.h>

namespace Lucee
{
/**
 * Base class for use in registering new creators.
 */
  template<class B>
  class ObjRegisterBase
  {
    public:
/**
 * Register a new object which will be created by its given name.
 *
 * @param nm Name by which object will be created.
 */
      ObjRegisterBase(const std::string& nm) 
        : name(nm)
      {
        Lucee::ObjCreator<B>::addCreator(name, this);
      }

/**
 * Delete the object, unregistering the creator class.
 */
      virtual ~ObjRegisterBase() 
      {
        Lucee::ObjCreator<B>::removeCreator(name);
      }

/**
 * Return a newly allocated object. Must be provided by derived
 * classes.
 *
 * @return Newly allocated object.
 */
      virtual B* getNew() = 0;

    private:
/** Name by which this object is to be created */
      std::string name;
  };

/**
 * Class to perform registration.
 */
  template<class D, class B>
  class ObjRegister : public ObjRegisterBase<B>
  {
    public:
/**
 * Register a new object which will be created by its given name.
 *
 * @param nm Name by which object will be created.
 */
      ObjRegister(const std::string& nm)
        : ObjRegisterBase<B>(nm)
      {
      }

/**
 * Return a newly allocated object. Must be provided by derived
 * classes.
 *
 * @return Newly allocated object.
 */
      B* getNew() 
      {
        D *d = new D;
        return dynamic_cast<B*>(d);
      }
  };

}

#endif // LC_OBJ_REGISTER_H
