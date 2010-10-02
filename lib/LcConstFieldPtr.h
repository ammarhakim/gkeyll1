/**
 * @file	LcConstFieldPtr.h
 *
 * @brief	Pointer to values stored in a field.
 *
 * @version	$Id: LcConstFieldPtr.h 222 2009-11-17 04:46:22Z a.hakim777 $
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

#ifndef LC_CONST_FIELD_PTR_H
#define LC_CONST_FIELD_PTR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <vector>

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
  class ConstFieldPtr
  {
    public:
/** Friend Field so it can create ptr objects */
      template <unsigned NDIM, typename TT> friend class Lucee::Field;

/**
 * Copy from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      ConstFieldPtr(const ConstFieldPtr<T>& ptr);

/**
 * Assign from a given field pointer.
 *
 * @param ptr Pointer to copy from.
 */
      ConstFieldPtr& operator=(const ConstFieldPtr<T>& ptr);

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
      const T& operator[](unsigned n) const { return data[n]; }

/**
 * Return pointer to underlying data.
 *
 * @return Pointer to data.
 */
      operator const T* () { return data; }

    private:
/**
 * Create a new field pointer to point to given data pointer.
 *
 * @param nc Number of components held by pointer.
 * @param data Data array. This should have length nc.
 */
      ConstFieldPtr(unsigned nc, const T *data);

/**
 * Set pointer to point to given data.
 *
 * @param dt Target data.
 */
      void setData(const T *dt) { data = dt; }

/** Number of components in field */
      unsigned numComponents;
/** Pointer to field data */
      const T *data;
  };
}

#endif  // LC_CONST_FIELD_PTR_H
