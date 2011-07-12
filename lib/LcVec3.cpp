/**
 * @file	LcVec3.cpp
 *
 * @brief	A vector in 3D space.
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
  template <typename T>
  Vec3<T>::Vec3(const T& val)
    : Lucee::FixedVector<3, T>(val) 
  {
  }

  template <typename T>
  Vec3<T>::Vec3(T x, T y, T z)
    : Lucee::FixedVector<3, T>(x, y, z)
  {
  }

  template <typename T>
  Vec3<T>::Vec3(const T xyz[3])
    : Lucee::FixedVector<3, T>(xyz)
  {
  }

  template <typename T>
  T
  Vec3<T>::dot(const Vec3<T>& vec) const
  {
    return vec[0]*this->getVal(0) + vec[1]*this->getVal(1) + vec[2]*this->getVal(2);
  }

  template <typename T>
  Vec3<T>
  Vec3<T>::cross(const Vec3<T>& vec) const
  {
    return Vec3<T>(
      this->getVal(1)*vec[2]-this->getVal(2)*vec[1],
      this->getVal(2)*vec[0]-this->getVal(0)*vec[2],
      this->getVal(0)*vec[1]-this->getVal(1)*vec[0]);
  }

  template <typename T>
  void
  Vec3<T>::normalize()
  {
    T len = ::sqrt(dot(*this));
    for (int i=0; i<3; ++i)
      setVal(i, this->getVal(i)/len);
  }

  template <typename T>
  Vec3<T>
  Vec3<T>::operator+(const Vec3<T>& vec) const
  {
    return Vec3<T>(
      vec[0]+this->getVal(0), vec[1]+this->getVal(1), vec[2]+this->getVal(2));
  }

  template <typename T>
  Vec3<T>
  Vec3<T>::operator-(const Vec3<T>& vec) const
  {
    return Vec3<T>(
      this->getVal(0)-vec[0], this->getVal(1)-vec[1], this->getVal(2)-vec[2]);
  }

  template <typename T>
  void
  Vec3<T>::scale(T fact)
  {
    for (unsigned i=0; i<3; ++i)
      setVal(i, this->getVal(i)*fact);
  }

// instantiations
  template class Vec3<float>;
  template class Vec3<double>;
}
