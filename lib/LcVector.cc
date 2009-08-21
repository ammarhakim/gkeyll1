/**
 * @file	LcVector.cc
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
        &Lucee::FixedVector<1,unsigned>(len)[0])
  {
  }

  template <typename T>
  Vector<T>::Vector(unsigned len, int start)
    : Lucee::Array<1, T>(
        &Lucee::FixedVector<1,unsigned>(len)[0], &Lucee::FixedVector<1,int>(start)[0])
  {
  }

  template <typename T>
  Vector<T>&
  Vector<T>::operator=(const T& val)
  {
    Lucee::Array<1, T>::operator=(val);
    return *this;
  }

// instantiations
  template class Lucee::Vector<int>;
  template class Lucee::Vector<float>;
  template class Lucee::Vector<double>;
}
