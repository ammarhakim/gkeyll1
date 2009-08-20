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
 * Get length of vector
 *
 * @return length of vector.
 */
      unsigned getLength() const { return this->getShape(0); }

  };
}

#endif // LC_VECTOR_H
