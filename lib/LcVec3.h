/**
 * @file	LcVec3.h
 *
 * @brief	A vector in 3D space.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
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
  class Vec3 : public Lucee::FixedVector<3, double>
  {
    public:
/**
 * Create new 3D vector with all components set to specified value
 *
 * @param val Value of all components. 
 */
      Vec3(const double& val=0.0);

/**
 * Create new 3D vector from specified values.
 *
 * @param x X-component of vector.
 * @param y Y-component of vector.
 * @param z Z-component of vector.
 */
      Vec3(double x, double y, double z);

/**
 * Construct new 3D vector from specified values.
 *
 * @param xyz Values of vector elements.
 */
      Vec3(const double xyz[3]);

/**
 * Compute the dot product of this vector with the supplied one.
 *
 * @param vec Vector to dot with.
 * @return dot product.
 */
      double dot(const Vec3& vec) const;

/**
 * Compute the cross product of this vector with the supplied
 * one. Computes res = this x vec.
 *
 * @param vec Vector to cross with.
 * @return cross product.
 */
      Vec3 cross(const Vec3& vec) const;

/**
 * Normalize the vector so its lenght it 1.0.
 */
      void normalize();

/**
 * Add this vector and supplied one and return new vector.
 *
 * @param vec Vector to add.
 * @return added vector.
 */
      Vec3 operator+(const Vec3& vec) const;

/**
 * Subtract supplied vector from this one.
 *
 * @param vec Vector to add.
 * @return subtracted vector.
 */
      Vec3 operator-(const Vec3& vec) const;

/**
 * Scale all elements with supplied factor.
 *
 * @param fact Factor to scale by.
 */
      void scale(double fact);
  };
}

#endif // LC_VEC_3_H
