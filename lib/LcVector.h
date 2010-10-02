/**
 * @file	LcVector.h
 *
 * @brief	Vector class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_VECTOR_H
#define LC_VECTOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>

namespace Lucee
{
/**
 * One-dimensional vector of data.
 */
  template <typename T>
  class Vector : public Lucee::Array<1, T>
  {
    public:
/**
 * Construct vector with specified length.
 *
 * @param len Length of the vector.
 */      
      Vector(unsigned len);

/**
 * Construct vector with specified length and start index.
 *
 * @param len Length of the vector.
 * @param start Start index for the vector.
 */
      Vector(unsigned len, int start);

/**
 * Create a new vector from input vector.
 *
 * @param vec Vector to create from.
 */
      Vector(const Vector<T>& vec);

/**
 * Create a new vector from 1D array object. Created vector and array
 * share data.
 *
 * @param arr Array to create from.
 */
      Vector(const Lucee::Array<1, T>& arr);

/**
 * Copy input vector.
 *
 * @param vec Vector to copy from.
 * @return Reference to this vector.
 */
      Vector<T>& operator=(const Vector<T>& vec);

/**
 * Assign all elements in vector to specified value.
 *
 * @param val Value to assign.
 * @return Reference to this vector.
 */
      Vector<T>& operator=(const T& val);

/**
 * Return value at i-th index.
 *
 * @param i index location.
 * @return value at location.
 */
      T operator[](int i) const { return this->operator()(i); }

/**
 * Return value at i-th index.
 *
 * @param i index location.
 * @return value at location.
 */
      T& operator[](int i) { return this->operator()(i); }

/**
 * Multiply vector by factor.
 *
 * @param fact Factor to multiply by.
 */
      void scale(const T& fact);

/**
 * Get length of vector
 *
 * @return length of vector.
 */
      unsigned getLength() const { return this->getShape(0); }

/**
 * Get length of vector
 *
 * @return length of vector.
 */
      unsigned size() const { return this->getShape(0); }

/**
 * Duplicate this vector.
 *
 * @return Copy of this vector.
 */
      Vector<T> duplicate() const;

/**
 * Return inner-product of vector with given vector.
 *
 * @param vec Vector to compute inner-product with.
 * @return inner-product.
 */
      T innerProduct(const Vector<T>& vec) const;

    private:
/**
 * Create a new vector of give size by reusing data pointer.
 *
 * @param len Length of vector.
 * @param dp Data pointer.
 */
      Vector(unsigned len, T *dp);
  };
}

#endif // LC_VECTOR_H
