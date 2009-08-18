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

namespace Lucee
{
  template <unsigned NELEM, typename T>
  class FixedVector
  {
    public:
/**
 * Construct new fixed-size vector setting all values to 0.
 */
      FixedVector(T init=T(0));

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

    private:
/** Data */
      T data[NELEM];
  };

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T init)
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = init;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T vals[NELEM])
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = vals[i];
  }
}

#endif // LC_FIXED_VECTOR_H
