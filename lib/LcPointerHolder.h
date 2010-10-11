/**
 * @file	LcPointerHolder.h
 *
 * @brief	Simple class to hold pointer
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_POINTER_HOLDER_H
#define LC_POINTER_HOLDER_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcExcept.h>
#include <LcLuaObjTypeId.h>

// lua includes
#include <lua.hpp>

// std includes
#include <typeinfo>

namespace Lucee
{
/** A simple class to wrap a pointer to an object */
  template <typename T>
  struct PointerHolder
  {
/**
 * Get a pointer to an object of type T stored as a Lua object. This
 * method is not always safe to use as there is no guarantee that the
 * cast from Lua stored void* to T is valid.
 *
 * @param L Lua state object to use for the user-data.
 * @return pointer to user object.
 */
      static T* getObj(lua_State *L)
      {
        Lucee::PointerHolder<T> *ph =
          (Lucee::PointerHolder<T>*) lua_touserdata(L, 1);
        return ph->pointer;
      }

/**
 * Get a pointer to an object of type T that is a base class for a
 * system of Lucee classes.
 *
 * @param L Lua state object to use for the user-data.
 * @param loc Location in Lua stack.
 * @return pointer to user object.
 */
      static T* getObjAsBase(lua_State *L, int loc=1)
      {
        void *obj;
        bool status = 
          Lucee::LuaObjTypeId::checkBaseTypeId(L, typeid(T).name(), &obj, loc);
        if (status == false)
        {
          Lucee::Except lce("Error fetching Lua userdata");
          throw lce;
        }
        PointerHolder<T> *ph = (PointerHolder<T>*) (obj);
        return ph->pointer;
      }

/**
 * Get a pointer to an object of type T that is a derived class of a
 * system of Lucee classes.
 *
 * @param L Lua state object to use for the user-data.
 * @param loc Location in Lua stack.
 * @return pointer to user object.
 */
      static T* getObjAsDerived(lua_State *L, int loc=1)
      {
        void *obj;
        bool status = 
          Lucee::LuaObjTypeId::checkDerivedTypeId(L, typeid(T).name(), &obj, loc);
        if (status == false)
        {
          Lucee::Except lce("Error fetching Lua userdata");
          throw lce;
        }
        PointerHolder<T> *ph = (PointerHolder<T>*) (obj);
        return ph->pointer;
      }

/**
 * Delete the held object.
 *
 * @param L Lua state object to use for the user-data.
 */
      static int deleteObject(lua_State *L)
      {
        T *myPtr = PointerHolder<T>::getObjAsDerived(L);
        delete myPtr;
        return 0;
      }

/** pointer to held object */
      T *pointer;
  };
}

#endif // LC_POINTER_HOLDER_H
