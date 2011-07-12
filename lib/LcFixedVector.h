/**
 * @file	LcFixedVector.h
 *
 * @brief	A fixed-sized Vector class.
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
/**
 * Vector of compiled time fixed-size data.
 */
  template <unsigned NELEM, typename T>
  class FixedVector
  {
    public:
/**
 * Construct new fixed-size vector with specified initial
 * values. Extra values are ignored and missing data is set to last
 * supplied value.
 *
 * @param v1 Value of first element.
 */
      FixedVector(T v1);

/**
 * Construct new fixed-size vector with specified initial
 * values. Extra values are ignored and missing data is set to last
 * supplied value.
 *
 * @param v1 Value of first element.
 * @param v2 Value of second element.
 */
      FixedVector(T v1, T v2);

/**
 * Construct new fixed-size vector with specified initial
 * values. Extra values are ignored and missing data is set to last
 * supplied value.
 *
 * @param v1 Value of first element.
 * @param v2 Value of second element.
 * @param v2 Value of thrid element.
 */
      FixedVector(T v1, T v2, T v3);

/**
 * Construct new fixed-size vector with specified initial values.
 *
 * @param vals Values of vector elements.
 */
      FixedVector(const T vals[NELEM]);

/**
 * Create new fixed-vector from existing one.
 *
 * @param fv Fixed-vector to copy from.
 */
      FixedVector(const FixedVector<NELEM, T>& fv);

/**
 * Create new fixed-vector from existing one.
 *
 * @param fv Fixed-vector to copy from.
 * @return reference to created vector.
 */
      FixedVector<NELEM, T>& operator=(const FixedVector<NELEM, T>& fv);

/**
 * Create new fixed-vector from an array.
 *
 * @param arr Array to copy from.
 * @return reference to created vector.
 */
      FixedVector<NELEM, T>& operator=(const T arr[NELEM]);

/**
 * Create new fixed-vector from a scalar value. Same value is assigned
 * to all elements.
 *
 * @param val Array to copy from.
 * @return reference to created vector.
 */
      FixedVector<NELEM, T>& operator=(const T& val);

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
      T getNorm() const;

    protected:
/**
 * Get value at location.
 *
 * @param i index location.
 * @return value.
 */
      const T& getVal(unsigned i) const { return data[i]; }

/**
 * Set value at location.
 *
 * @param i index location.
 * @param val value to set
 */
      void setVal(unsigned i, const T& val) { data[i] = val; }

    private:
/** Data */
      T data[NELEM];
  };

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T v1)
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = v1;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T v1, T v2)
  {
      switch (NELEM)
      {
          case 0:
            break;
          case 1:
            data[0] = v1;
            break;
          default:
            data[0] = v1;
            data[1] = v2;
            break;
      }
      for (unsigned i=2; i<NELEM; ++i)
        data[i] = v2;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(T v1, T v2, T v3)
  {
      switch (NELEM)
      {
          case 0:
            break;
          case 1:
            data[0] = v1;
            break;
          case 2:
            data[0] = v1;
            data[1] = v2;
            break;
          default:
            data[0] = v1;
            data[1] = v2;
            data[2] = v3;
            break;
      }
      for (unsigned i=3; i<NELEM; ++i)
        data[i] = v2;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(const T vals[NELEM])
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = vals[i];
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>::FixedVector(const FixedVector<NELEM, T>& fv)
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = fv.data[i];
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>& 
  FixedVector<NELEM, T>::operator=(const FixedVector<NELEM, T>& fv)
  {
    if (&fv == this)
      return *this;

    for (unsigned i=0; i<NELEM; ++i)
      data[i] = fv.data[i];
    return *this;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>& 
  FixedVector<NELEM, T>::operator=(const T arr[NELEM])
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = arr[i];
    return *this;
  }

  template <unsigned NELEM, typename T>
  FixedVector<NELEM, T>&
  FixedVector<NELEM, T>::operator=(const T& val)
  {
    for (unsigned i=0; i<NELEM; ++i)
      data[i] = val;
    return *this;
  }

  template <unsigned NELEM, typename T>
  T
  FixedVector<NELEM, T>::getNorm() const
  {
    T norm = 0;
    for (unsigned i=0; i<NELEM; ++i)
      norm += data[i]*data[i];
    return std::sqrt(norm);
  }
}

#endif // LC_FIXED_VECTOR_H
