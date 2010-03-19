/**
 * @file	LcMatrix.cpp
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
/** Transpose mask */
#define LC_TRANSPOSE matrixMasks[0]
/** Set transpose mask */
#define LC_SET_TRANSPOSE(bit) (bit) |= LC_TRANSPOSE
/** Clear transpose mask */
#define LC_CLEAR_TRANSPOSE(bit) (bit) &= ~LC_TRANSPOSE
/** Check if transpose mask is set */
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
  Lucee::Vector<T>
  Matrix<T>::getRow(unsigned row)
  {
    if (row > numRows()-1)
    {
      Lucee::Except lce("Matrix::getRow: Incorrect row ");
      lce << row << " specified"  << std::endl;
      throw lce;
    }

    unsigned defDims[1] = {0};
    int defDimsIdx[1];
    defDimsIdx[0] = row;
// call parent class to get row
    Lucee::Array<1, T> myRow = this->template deflate<1>(defDims, defDimsIdx);
    return Lucee::Vector<T>(myRow);
  }

  template <typename T>
  Lucee::Vector<T>
  Matrix<T>::getCol(unsigned col)
  {
    if (col > numColumns()-1)
    {
      Lucee::Except lce("Matrix::getCol: Incorrect col ");
      lce << col << " specified"  << std::endl;
      throw lce;
    }

    unsigned defDims[1] = {1};
    int defDimsIdx[1];
    defDimsIdx[0] = col;
// call parent class to get column
    Lucee::Array<1, T> myCol = this->template deflate<1>(defDims, defDimsIdx);
    return Lucee::Vector<T>(myCol);
  }

  template <typename T>
  void
  Matrix<T>::scaleRows(const Lucee::Vector<double>& fv)
  {
    if (fv.getLength() != numRows())
      throw Lucee::Except("Matrix::scaleRows: Number of elements should match number of rows");

    int ifv = fv.getLower(0);
    for (int i=this->getLower(0); i<this->getUpper(0); ++i)
    {
      double t1 = fv[ifv++];
      for (int j=this->getLower(1); j<this->getUpper(1); ++j)
        this->operator()(i,j) *= t1;
    }
  }

  template <typename T>
  void
  Matrix<T>::scaleCols(const Lucee::Vector<double>& fv)
  {
    if (fv.getLength() != numColumns())
      throw Lucee::Except("Matrix::scaleCols: Number of elements should match number of columns");

    int ifv = fv.getLower(0);
    for (int j=this->getLower(1); j<this->getUpper(1); ++j)
    {
      double t1 = fv[ifv++];
      for (int i=this->getLower(0); i<this->getUpper(0); ++i)
        this->operator()(i,j) *= t1;
    }
  }

  template <typename T>
  Lucee::Matrix<T>
  Matrix<T>::getView(int lr, int ur, int lc, int uc)
  {
    int lo[2], up[2];
    lo[0] = lr; lo[1] = lc;
    up[0] = ur; up[1] = uc;
    Lucee::Region<2, int> rgn(lo, up);
    Lucee::Array<2, T, Lucee::ColMajorIndexer> slice
      = this->getSlice(rgn);
    return Matrix<T>(slice);
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

  template <typename T>
  Matrix<T>::Matrix(Lucee::Array<2, T, Lucee::ColMajorIndexer>& arr)
    : Lucee::Array<2, T, Lucee::ColMajorIndexer>(arr)
  {
  }

// instantiations
  template class Lucee::Matrix<int>;
  template class Lucee::Matrix<float>;
  template class Lucee::Matrix<double>;
}
