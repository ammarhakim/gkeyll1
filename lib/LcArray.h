/**
 * @file	LcArray.h
 *
 * @brief	Serial array class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_ARRAY_H
#define LC_ARRAY_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcObject.h>

namespace Lucee
{
  template <unsigned NDIM, typename T>
  class Array : Lucee::Object
  {
    public:
/**
 * Construct array with specified shape.
 *
 * @param shape Shape of the array.
 * @param init Initial value to assign to all elements.
 */      
      Array(unsigned shape[NDIM], const T& init=(T)0);

/**
 * Construct array with specified shape and start indices.
 *
 * @param shape Shape of the array.
 * @param start indices for the array.
 * @param init Initial value to assign to all elements.
 */
      Array(unsigned shape[NDIM], int start[NDIM], const T& init=(T)0);

/**
 * Get shape of array.
 *
 * @param shape On return contains the shape of the array.
 */
      void getShape(unsigned shape[NDIM]) const;

/**
 * Get shape of array in the specified direction.
 *
 * @return shape in specified direction.
 */
      unsigned getShape(unsigned dir) const;

/**
 * Get start index in the specified direction.
 *
 * @return start index in specified direction.
 */
      int getStart(unsigned dir) const;

/**
 * Get one past the end index in the specified direction.
 *
 * @return one past the end index in specified direction.
 */
      int getEnd(unsigned dir) const;

    private:
/** Shape of array */
      unsigned shape[NDIM];
/** Start index of array */
      int start[NDIM];
  };
}

#endif // LC_ARRAY_H
