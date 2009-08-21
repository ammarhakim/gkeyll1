/**
 * @file	LcFixedVector.h
 *
 * @brief	A fixed-sized Vector class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_FIXED_VECTOR_H
#define LC_FIXED_VECTOR_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// std includes
#include <stdarg.h>
#include <cmath>

namespace Lucee
{
  template <unsigned NELEM, typename T>
  class FixedVector
  {
    public:
/**
 * Construct new fixed-size vector with specified initial values.
 *
 * @param init Initial value to apply to all elements.
 */
      FixedVector(const T& init);

/**
 * Construct new fixed-size vector with specified initial values.
 *
 * @param vals v1 Value of first element.
 * @param vals v2 Value of second element.
 */
      FixedVector(T v1, T v2, ...);

/**
 * Construct new fixed-size vector with specified initial values.
 *
 * @param vals Values of vector elements.
 */
      FixedVector(T vals[NELEM]);

/**
 * Return number of elements in vector.
 *
 * @return number of elements in vector.
 */
      unsigned numElem() const { return NELEM; }

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

/**
 * Return 2-norm of vector.
 *
 * @return 2-norm of vector.
 */
      T get2Norm() const;

    private:
/** Data */
      T data[NELEM];
  };

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(const T& init)
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = init;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T v1, T v2, ...)
  {
    va_list elems;
    va_start(elems, v2);
    data[0] = v1;
    data[1] = v2;
    for (unsigned i=2; i<NELEM; ++i)
      data[i] = va_arg(elems, T);

    va_end(elems);
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T vals[NELEM])
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = vals[i];
  }

  template <unsigned NELEM, typename T>
  T
  FixedVector<NELEM, T>::get2Norm() const
  {
    T norm = 0.0;
    for (unsigned i=0; i<NELEM; ++i)
      norm += data[i]*data[i];
    return sqrt(norm);
  }
}

#endif // LC_FIXED_VECTOR_H
