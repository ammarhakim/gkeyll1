/**
 * @file	LcVector.cpp
 *
 * @brief	Vector class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcVector.h>
#include <LcFixedVector.h>

namespace Lucee
{
  template <typename T>
  Vector<T>::Vector(unsigned len)
    : Lucee::Array<1, T>(
      &Lucee::FixedVector<1, unsigned>(len)[0])
  {
  }

  template <typename T>
  Vector<T>::Vector(unsigned len, int start)
    : Lucee::Array<1, T>(
      &Lucee::FixedVector<1, unsigned>(len)[0], &Lucee::FixedVector<1,int>(start)[0])
  {
  }

  template <typename T>
  Vector<T>::Vector(const Vector<T>& vec)
    : Lucee::Array<1, T>(vec)
  {
  }

  template <typename T>
  Vector<T>::Vector(const Lucee::Array<1, T>& arr)
    : Lucee::Array<1, T>(arr)
  {
  }

  template <typename T>
  Vector<T>&
  Vector<T>::operator=(const Vector<T>& vec)
  {
    if (&vec == this)
      return *this;
    Lucee::Array<1, T>::operator=(vec);
    return *this;
  }

  template <typename T>
  Vector<T>&
  Vector<T>::operator=(const T& val)
  {
    Lucee::Array<1, T>::operator=(val);
    return *this;
  }

  template <typename T>
  void
  Vector<T>::scale(const T& fact)
  {
    for (int i=this->getLower(0); i<this->getUpper(0); ++i)
      this->operator[](i) *= fact;
  }

  template <typename T>
  Vector<T>
  Vector<T>::duplicate() const
  {
    unsigned shape[1];
    int start[1];
// get out shape and start indices
    this->fillWithShape(shape);
    this->fillWithStart(start);
// create new vector and copy data into it
    Vector<T> dup(shape[0], start[0]);
    for (int i=this->getLower(0); i<this->getUpper(0); ++i)
      dup[i] = this->operator[](i);
    return dup;
  }

  template <typename T>
  Vector<T>::Vector(unsigned len, T *dp)
    : Lucee::Array<1, T>(
      Lucee::Region<1, int>(&Lucee::FixedVector<1, int>(len)[0]),
      dp)
  {
  }

// instantiations
  template class Lucee::Vector<int>;
  template class Lucee::Vector<float>;
  template class Lucee::Vector<double>;
}
