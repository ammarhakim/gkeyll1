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

// lua includes
#include <lua.hpp>

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

namespace Lucee
{
/** A simple class to wrap a pointer to an object */
  template <typename T>
  struct PointerHolder
  {
/**
 * Check user type and if correct return pointer to help object.
 *
 * @param L Lua state object to use for the user-data.
 */
      static T* checkUserType(lua_State *L)
      {
        Lucee::PointerHolder<T> *ph =
          (Lucee::PointerHolder<T>*) luaL_checkudata(L, 1, typeid(T).name());
        return ph->pointer;
      }

/**
 * Delete the held object.
 *
 * @param L Lua state object to use for the user-data.
 */
      static int deleteObject(lua_State *L)
      {
        T *myPtr = PointerHolder<T>::checkUserType(L);
        delete myPtr;
        return 0;
      }

/** pointer to held object */
      T *pointer;
  };
}

#endif // LC_POINTER_HOLDER_H
