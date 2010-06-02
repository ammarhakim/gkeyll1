/**
 * @file	LcVec3.cpp
 *
 * @brief	A vector in 3D space.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2010, Ammar Hakim.
 */

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcVec3.h>

// std includes
#include <cmath>

namespace Lucee
{
  Vec3::Vec3(double x, double y, double z)
    : Lucee::FixedVector<3, double>(x, y, z) 
  {
  }

  Vec3::Vec3(const double xyz[3])
    : Lucee::FixedVector<3, double>(xyz)
  {
  }

  double
  Vec3::dot(const Vec3& vec) const
  {
    return vec[0]*getVal(0) + vec[1]*getVal(1) + vec[2]*getVal(2);
  }

  Vec3
  Vec3::cross(const Vec3& vec) const
  {
    return Vec3(
      getVal(1)*vec[2]-getVal(2)*vec[1],
      getVal(2)*vec[0]-getVal(0)*vec[2],
      getVal(0)*vec[1]-getVal(1)*vec[0]);
  }

  void
  Vec3::normalize()
  {
    double len = ::sqrt(dot(*this));
    for (int i=0; i<3; ++i)
      setVal(i, getVal(i)/len);
  }

  Vec3
  Vec3::operator+(const Vec3& vec) const
  {
    return Vec3(
      vec[0]+getVal(0), vec[1]+getVal(1), vec[2]+getVal(2));
  }

  Vec3
  Vec3::operator-(const Vec3& vec) const
  {
    return Vec3(
      getVal(0)-vec[0], getVal(1)-vec[1], getVal(2)-vec[2]);
  }

  void
  Vec3::scale(double fact)
  {
    for (unsigned i=0; i<3; ++i)
      setVal(i, getVal(i)*fact);
  }
}
