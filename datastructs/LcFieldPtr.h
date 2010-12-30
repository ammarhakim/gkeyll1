/**
 * @file	LcFieldPtr.h
 *
 * @brief	Pointer to values stored in a field.
 *
 * @version	$Id: LcFieldPtr.h 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_FIELD_PTR_H
#define LC_FIELD_PTR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcVector.h>

// std includes
#include <vector>

namespace Lucee
{
// forward declaration for making Field class friend
  template <unsigned NDIM, typename TT> class Field;
// forward declaration for making ConstFieldPtr class friend
  template <typename TT> class ConstFieldPtr;

/**
 * A field pointer can be used to access/modify the elements of a
 * field. This is a private class that can only be created by
 * Lucee:Field objects.
 */
  template <typename T>
  class FieldPtr
  {
    public:
/** Friend Field so it can create ptr objects */
      template <unsigned NDIM, typename TT> friend class Lucee::Field;

/** Friend ConstFieldPtr so it can create ptr objects */
      template <typename TT> friend class Lucee::ConstFieldPtr;

/**
 * Copy from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      FieldPtr(const FieldPtr<T>& ptr);

/**
 * Create from a std::vector. The new object shares data with the
 * supplied vector. Hence, the created object is no longer valid when
 * the supplied std::vector goes away.
 *
 * @param vec Vector to create from.
 */
      FieldPtr(std::vector<T>& vec);

/**
 * Create from a Lucee::Vector object. The new object shares data with
 * the supplied vector. Hence, the created object is no longer valid when
 * the supplied vector goes away.
 *
 * @param vec Vector to create from.
 */
      FieldPtr(Lucee::Vector<T> vec);

/**
 * Create an empty field pointer that allows storing data.
 *
 * @param num Number of elements in vector.
 */
      FieldPtr(unsigned num);

/**
 * Delete field pointer.
 */
      ~FieldPtr();

/**
 * Assign from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      FieldPtr& operator=(const FieldPtr<T>& ptr);

/**
 * Assign all values to supplied one.
 *
 * @param val Value to set.
 */
      FieldPtr& operator=(const T& val);

/**
 * Number of elements indexed by pointer.
 *
 * @return Number of elements indexed by pointer.
 */
      unsigned getNumComponents() const { return numComponents; }

/**
 * Return element in field data.
 *
 * @param n Element index to fetch.
 */
      T& operator[](int n) { return data[n]; }

/**
 * Return element in field data.
 *
 * @param n Element index to fetch.
 */
      T operator[](int n) const { return data[n]; }

/**
 * Return pointer to underlying data.
 *
 * @return Pointer to data.
 */
      operator T* () { return data; }

    private:
/**
 * Create a new field pointer to point to given data pointer.
 *
 * @param nc Number of components held by pointer.
 * @param data Data array. This should have length nc.
 */
      FieldPtr(unsigned nc, T *data);

/**
 * Set pointer to point to given data.
 *
 * @param dt Target data.
 */
      void setData(T *dt) { data = dt; }

/** Number of components in field */
      unsigned numComponents;
/** Pointer to field data */
      T *data;
/** Was data allocated */
      bool isAlloc;
  };
}

#endif  // LC_FIELD_PTR_H
