/**
 * @file	LcRegisteredObjList.h
 *
 * @brief	Class to hold list of registered objects.
 */

#ifndef LC_REGISTERED_OBJ_LIST_H
#define LC_REGISTERED_OBJ_LIST_H

// lucee includes
#include <LcObjRegistryBase.h>
#include <LcObjRegistry.h>

// std includes
#include <vector>

namespace Lucee
{
/**
 * This class holds list of objects registered for a specified base
 * class. Allows cleaning up of memory in a systematic way.
 */
  template <typename B>
  class RegisteredObjList
  {
    public:
/**
 * Delete all registered objects.
 */
      ~RegisteredObjList()
      {
        std::vector<Lucee::ObjRegistryBase<B>*>::iterator itr
          = objList.begin();
        for ( ; itr != objList.end(); ++itr)
          delete *itr;
        objList.clear();
      }

/**
 * Append a new object to list.
 *
 * @return reference to this object.
 */
      template <typename D>
      RegisteredObjList<B>& append()
      {
        objList.push_back( new ObjRegistry<B,D> );
        return *this;
      }

    private:
/** List of objects registered with base class B */
      std::vector<Lucee::ObjRegistryBase<B>*> objList;
  }
}

#endif // LC_REGISTERED_OBJ_LIST_H
