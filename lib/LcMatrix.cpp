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
    : Lucee::Array<2, T, Lucee::ColMajorIndexer>(
        &Lucee::FixedVector<2, unsigned>(row, col)[0])
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer>(shape)
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2], int start[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer>(shape, start)
  {
    LC_CLEAR_TRANSPOSE(traits);
  }

  template <typename T>
  Matrix<T>::Matrix(const Matrix<T>& mat)
    : Lucee::Array<2, T, Lucee::ColMajorIndexer>(mat)
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
    Lucee::Array<2, T, Lucee::ColMajorIndexer>::operator=(mat);
    return *this;
  }

  template <typename T>
  Matrix<T>&
  Matrix<T>::operator=(const T& val)
  {
    Lucee::Array<2, T, Lucee::ColMajorIndexer>::operator=(val);
    return *this;
  }

  template <typename T>
  Matrix<T>
  Matrix<T>::duplicate() const
  {
    unsigned shape[2];
    int start[2];
// get our shape and start indices
    this->fillWithShape(shape);
    this->fillWithStart(start);
// create new matrix and copy data into it
    Matrix<T> dup(shape, start);
    for (int i=this->getLower(0); i<this->getUpper(0); ++i)
      for (int j=this->getLower(1); j<this->getUpper(1); ++j)
        dup(i,j) = this->operator()(i,j);
    return dup;
  }

  template <typename T>
  bool
  Matrix<T>::isTranspose() const
  {
    return LC_IS_TRANSPOSE(traits);
  }

// instantiations
  template class Lucee::Matrix<int>;
  template class Lucee::Matrix<float>;
  template class Lucee::Matrix<double>;
}
