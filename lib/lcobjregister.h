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

#ifndef LC_OBJ_REGISTER
#define LC_OBJ_REGISTER

// lib includes
#include <lcexcept.h>
#include <lcobjcreator.h>

namespace Lucee
{

/**
 * 
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
      ObjRegister(const std::string& name)
        : ObjRegisterBase<B>(name) 
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

#endif // LC_OBJ_REGISTER
