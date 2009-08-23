/**
 * @file	LcMatrix.cc
 *
 * @brief	Matrix class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcFixedVector.h>
#include <LcMatrix.h>

// std include
#include <vector>

namespace Lucee
{
// Masks for matrix traits
  static unsigned int matrixMasks[] =
  {0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80};

// set of macros for setting/getting matrix traits
#define LC_TRANSPOSE matrixMasks[0]
#define LC_SET_TRANSPOSE(bit) (bit) |= LC_TRANSPOSE
#define LC_CLEAR_TRANSPOSE(bit) (bit) &= ~LC_TRANSPOSE
#define LC_IS_TRANSPOSE(bit) (bit) & LC_TRANSPOSE

  template <typename T>
  Matrix<T>::Matrix(unsigned row, unsigned col)
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(
        &Lucee::FixedVector<2, unsigned>(row, col)[0])
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(shape)
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2], int start[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(shape, start)
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(const Matrix<T>& mat)
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(mat)
  {
    traits = mat.traits;
  }

  template <typename T>
  Matrix<T>&
  Matrix<T>::operator=(const Matrix<T>& mat)
  {
    if (&mat == this)
      return *this;

    traits = mat.traits;
// call base class assignment operator
    Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >::operator=(mat);
    return *this;
  }

  template <typename T>
  Matrix<T>&
  Matrix<T>::operator=(const T& val)
  {
    Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >::operator=(val);
    return *this;
  }

// instantiations
  template class Lucee::Matrix<int>;
  template class Lucee::Matrix<float>;
  template class Lucee::Matrix<double>;
}
