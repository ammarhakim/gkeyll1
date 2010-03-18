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

namespace Lucee
{
// forward declaration for making Field class friend
  template <unsigned NDIM, typename TT> class Field;

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

/**
 * Copy from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      FieldPtr(const FieldPtr<T>& ptr);

/**
 * Assign from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      FieldPtr& operator=(const FieldPtr<T>& ptr);

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
      T& operator[](unsigned n) { return data[n]; }

/**
 * Return element in field data.
 *
 * @param n Element index to fetch.
 */
      T operator[](unsigned n) const { return data[n]; }

/**
 * Get a Lucee::Vector class from this field pointer.
 *
 * @param vector class for pointer.
 */
      Lucee::Vector<T> asVector();

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
  };
}

#endif  // LC_FIELD_PTR_H
