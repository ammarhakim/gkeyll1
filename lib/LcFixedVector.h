/**
 * @file	LcFixedVector.h
 *
 * @brief	A fixed-sized Vector class.
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

namespace Lucee
{
  template <usigned NELEM, typename T>
  class FixedVector
  {
    public:
/**
 * Construct new fixed-size vector setting all values to 0.
 */
      FixedVector(const T& init=T(0));

/**
 * Construct new fixed-size vector with specified initial values.
 *
 * @param vals Values of vector elements.
 */
      FixedVector(T vals[NELEM]);

/**
 * Return value of vector at location.
 *
 * @param i index location.
 * @return value.
 */
      const T& operator[](unsigned i) const { return data[i]; }

/**
 * Return settable value of vector at location.
 *
 * @param i index location.
 * @return value.
 */
      T& operator[](unsigned i) { return data[i]; }

    private:
/** Data stored in fixed-size vector */
      T data[NELEM];
  };
}

#endif // LC_VECTOR_H
