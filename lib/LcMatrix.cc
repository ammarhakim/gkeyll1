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
#include <LcFixedVector.h>
#include <LcMatrix.h>

namespace Lucee
{
  template <typename T>
  Matrix<T>::Matrix(unsigned row, unsigned col)
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(
        &Lucee::FixedVector<2, unsigned>(row, col)[0])
  {
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(shape)
  {
  }

  template <typename T>
  Matrix<T>::Matrix(unsigned shape[2], int start[2])
    : Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >(shape, start)
  {
  }

  template <>
  void
  Matrix<double>::eig(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi)
  {
  }

  template <>
  void 
  Matrix<double>::eig(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi, 
    Matrix<double>& vecl, Matrix<double>& vecr)
  {
  }

  template <>
  void
  Matrix<double>::eigRight(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi, 
    Matrix<double>& vec)
  {
  }

  template <>
  void 
  Matrix<double>::eigLeft(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi,
    Matrix<double>& vec)
  {
  }

  template <>
  void
  Matrix<double>::solve(Lucee::Matrix<double>& rhs)
  {
  }

// instantiations
  template class Lucee::Matrix<int>;
  template class Lucee::Matrix<float>;
  template class Lucee::Matrix<double>;
}
