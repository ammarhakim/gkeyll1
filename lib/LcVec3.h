/**
 * @file	LcVec3.h
 *
 * @brief	A vector in 3D space.
 */

#ifndef LC_VEC_3_H
#define LC_VEC_3_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcFixedVector.h>

namespace Lucee
{
/**
 * Vector of doubles in 3D space. Provides special methods for dot,
 * cross products, creation of unit vectors etc.
 */
  template <typename T>
  class Vec3 : public Lucee::FixedVector<3, T>
  {
    public:
/**
 * Create new 3D vector with all components set to specified value
 *
 * @param val Value of all components. 
 */
      Vec3(const T& val=0.0);

/**
 * Create new 3D vector from specified values.
 *
 * @param x X-component of vector.
 * @param y Y-component of vector.
 * @param z Z-component of vector.
 */
      Vec3(T x, T y, T z);

/**
 * Construct new 3D vector from specified values.
 *
 * @param xyz Values of vector elements.
 */
      Vec3(const T xyz[3]);

/**
 * Compute the dot product of this vector with the supplied one.
 *
 * @param vec Vector to dot with.
 * @return dot product.
 */
      T dot(const Vec3<T>& vec) const;

/**
 * Compute the cross product of this vector with the supplied
 * one. Computes res = this x vec.
 *
 * @param vec Vector to cross with.
 * @return cross product.
 */
      Vec3<T> cross(const Vec3<T>& vec) const;

/**
 * Normalize the vector so its length it 1.0.
 */
      void normalize();

/**
 * Add this vector and supplied one and return new vector.
 *
 * @param vec Vector to add.
 * @return added vector.
 */
      Vec3<T> operator+(const Vec3<T>& vec) const;

/**
 * Subtract supplied vector from this one.
 *
 * @param vec Vector to add.
 * @return subtracted vector.
 */
      Vec3<T> operator-(const Vec3<T>& vec) const;

/**
 * Scale all elements with supplied factor.
 *
 * @param fact Factor to scale by.
 */
      void scale(T fact);
  };
}

#endif // LC_VEC_3_H
