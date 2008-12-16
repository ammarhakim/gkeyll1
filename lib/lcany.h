/**
 * @file	lcany.h
 *
 * @brief	Class to hold objects of arbitary type
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim. Released under Eclipse
 * Licence version 1.0.
 */

#ifndef LC_ANY_H
#define LC_ANY_H

// std includes
#include <algorithm>
#include <typeinfo>

namespace Lucee
{
/**
 * Class Any is based on the "any" class described in "Valued
 * Conversion", Kevlin Henney, C++ Report, July-August 2000, pages
 * 37--40. Any should be used for objects with simple types which
 * support a copy ctor.
 *
 * Example usage is:
 *
 * Any any(22); // store integer 22
 * cout << any_cast<int>(any) << endl; // prints 22
 *
 *  any = Any(string("Hello World"); // store string
 *  cout << any_cast<string>(any) << endl; // prints "Hello World"
 */
  class Any
  {
    public:
/**
 * Create empty object
 */
      Any()
        : content(0) 
      {
      }

/**
 * Create new Any from an object of a given type
 *
 * @param value Value of object stored.
 */
      template<typename VALUETYPE>
      Any(const VALUETYPE& value)
        : content(new _Holder<VALUETYPE>(value)) 
      {
      }

/**
 * Copy constructor: a deep copy is made
 *
 * @param other Any object to copy
 */
      Any(const Any& other)
        : content(other.content ? other.content->clone() : 0) 
      {
      }

/**
 * Destroy object
 */
      ~Any() 
      {
        // destroy held object
        delete content;
      }

/**
 * Swap contents of object with contents of supplied object
 *
 * @param rhs Replace the value in this Any with 'rhs'
 * @return Reference to this object
 */
      Any& swap(Any& rhs) 
      {
        std::swap(content, rhs.content);
        return *this;
      }

/**
 * Assignment operator: use Any object to create a new object
 *
 * @param rhs Any to assign from
 * @return Reference to this object
 */
      Any& operator=(const Any& rhs) 
      {
        Any(rhs).swap(*this);
        return *this;
      }

/**
 * Assignment operator: use VALUETYPE object to create a new object
 * 
 * @return Reference to this object
 */
      template<typename VALUETYPE>
      Any& operator=(const VALUETYPE& rhs) 
      {
        Any(rhs).swap(*this);
        return *this;
      }

/**
 * Is object empty?
 *
 * @return true if object is empty, false otherwise
 */
      bool isEmpty() const 
      {
        return !content;
      }

/**
 * Type_info for this object's held value
 *
 * @return type_info for held object
 */
      const std::type_info& type() const 
      {
        return content ? content->type() : typeid(void);
      }

/**
 * Convert the held object to a pointer. Do not use this function
 * directly without checking if the returned pointer is not NULL.
 *
 * @return pointer to held object
 */
      template<typename VALUETYPE>
      const VALUETYPE* to_ptr() const 
      {
        // return pointer to held object or NULL if wrong type
        return 
          type() == typeid(VALUETYPE) 
          ? &static_cast<_Holder<VALUETYPE> *> (content)->held
          : 0
          ;
      }

/**
 * Convert help object to a void *. This is not a safe operation as
 * it breaks typechecking.
 *
 * @return pointer to held object, with type information removed
 */
      const void* to_void_ptr() const 
      {
        return content->void_ptr();
      }

/**
 * Extract the data from the Any object. 
 *
 * If the type specified by the template 'VALUETYPE' is not the
 * correct type of the object stored a std::bad_cast exception is
 * thrown.
 *
 * @return value of held object
 */
      template<typename VALUETYPE>
      VALUETYPE to_value() const 
      {
        // return object or throw exception if types do not match
        const VALUETYPE * result =
          type() == typeid(VALUETYPE) 
          ? &static_cast<_Holder<VALUETYPE> *> (content)->held
          : 0
          ;
        // if it is not null, return value else thow exception
        if (result)
          return *result;
        else
          throw std::bad_cast();  
      }

    private:

/**
 * Internal base-class to hold the object
 */
      class _PlaceHolder
      {
        public:
/**
 * Destructor
 */
          virtual ~_PlaceHolder() 
          {
          }

/**
 * Type_info for this object's held value
 *
 * @return type_info for held object
 */
          virtual const std::type_info& type() const = 0;

/**
 * Make a copy: must be provided by children
 *
 * @return copy of wrapped object
 */
          virtual _PlaceHolder* clone() const = 0;

/**
 * Convert help object to a void *.
 *
 * @return pointer to held object, with type information removed
 */
          virtual const void * void_ptr() const = 0;
      };

/**
 * Internal class to hold the object
 */
      template<typename VALUETYPE>
      class _Holder : public _PlaceHolder
      {
        public:
/**
 * Create new object from supplied value
 *
 * @param value Value of object
 */
          _Holder(const VALUETYPE& value)
            : held(value) 
          {
          }

/**
 * Type_info for this object's held value
 *
 * @return type_info for held object
 */
          virtual const std::type_info& type() const 
          {
            return typeid(VALUETYPE);
          }

/**
 * Make a copy: must be provided by children
 *
 * @return Pointer to wrapped held object
 */
          virtual _PlaceHolder* clone() const 
          {
            return new _Holder(held);
          }

/**
 * Convert help object to a void *.
 *
 * @return pointer to held object, with type information removed
 */
          virtual const void * void_ptr() const 
          {
            return (const void*) &held;
          }

/** Value of held object */
          VALUETYPE held;
      };

/** Pointer to held object */
      _PlaceHolder *content;
  };

/**
 * Extract the data from the Any object. 
 *
 * @param operand Any object from which to extract value. If the
 * type specified by the template 'VALUETYPE' is not the correct type
 * of the object stored a std::bad_cast exception is thrown.
 */
  template<typename VALUETYPE>
  VALUETYPE
  any_cast(const Any& operand)
  {
    // return object or throw exception if types do not match
    const VALUETYPE * result = operand.template to_ptr<VALUETYPE>();

    // if it is not null, return value else thow exception
    if (result)
      return *result;
    else
      throw std::bad_cast();
  }

/**
 * Extract the data from the Any object. 
 *
 * @param operand Any object from which to extract value. If the
 * type specified by the template 'VALUETYPE' is not the correct type
 * of the object stored a std::bad_cast exception is thrown.
 */
  template<typename VALUETYPE>
  VALUETYPE
  any_cast(Any& operand)
  {
    // return object or throw exception if types do not match
    const VALUETYPE * result = operand.template to_ptr<VALUETYPE>();

    // if it is not null, return value else thow exception
    if (result)
      return *result;
    else
      throw std::bad_cast();
  }
}

#endif // LC_ANY_H
