/**
 * @file	LcObjRegistryBase.h
 *
 * @brief	Base class for registration system.
 */

#ifndef LC_OBJ_REGISTRY_BASE_H
#define LC_OBJ_REGISTRY_BASE_H

namespace Lucee
{
/**
 * This class provides a base class for the ObjRegistry class. It is
 * needed to help in managing registered object so memory leaks are
 * avoided.
 */
  template <class B>
  class ObjRegistryBase
  {
    public:
/**
 * Destroy object.
 */
      virtual ~ObjRegistryBase()
      {
// nothing to do
      }
  };
}

#endif // LC_OBJ_REGISTRY_BASE_H
