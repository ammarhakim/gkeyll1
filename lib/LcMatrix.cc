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
#include <LcLapackDeclarations.h>
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

  template <>
  void
  Matrix<double>::eig(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi)
  {
// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(numRows()*numColumns());
    unsigned count = 0;
    for (int j=getLower(1); j<getUpper(1); ++j)
      for (int i=getLower(0); i<getUpper(0); ++i)
        A[count++] = this->operator()(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = 1;
    LDVR = 1;

    LDA = numRows();
    N = LDA;
    JOBVR = 'N';
    JOBVL = 'N';
    LWORK = 5*N;
    std::vector<double> WORK(LWORK);

// call LAPACK routine to compute eigenvalues
    dgeev_(&JOBVL, &JOBVR, &N, &A[0], &LDA,
      &evr[evr.getLower(0)],
      &evi[evi.getLower(0)],
      0, &LDVL,
      0, &LDVR,
      &WORK[0], &LWORK,
      &INFO);

// check if eigensolver worked
    if (INFO!=0)
      throw Lucee::Except("Matrix::eig: Eigenvalue solver failed");
  }

  template <>
  void 
  Matrix<double>::eig(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi, 
    Matrix<double>& vecl, Matrix<double>& vecr)
  {
// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(numRows()*numColumns());
    unsigned count = 0;
    for (int j=getLower(1); j<getUpper(1); ++j)
      for (int i=getLower(0); i<getUpper(0); ++i)
        A[count++] = this->operator()(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = vecl.numRows();
    LDVR = vecl.numRows();

    LDA = numRows();
    N = LDA;
    JOBVR = 'V';
    JOBVL = 'V';
    LWORK = 5*N;
    std::vector<double> WORK(LWORK);
     
// call LAPACK routine to compute eigenvalues and eigenvectors
    dgeev_(&JOBVL, &JOBVR, &N, &A[0], &LDA,
      &evr[evr.getLower(0)],
      &evi[evi.getLower(0)],
      &vecl(vecl.getLower(0), vecl.getLower(1)), &LDVL,
      &vecr(vecr.getLower(0), vecr.getLower(1)), &LDVR,
      &WORK[0], &LWORK,
      &INFO);

// check if eigensolver worked
    if (INFO!=0)
      throw Lucee::Except("Matrix::eig: Eigenvalue solver failed");
  }

  template <>
  void
  Matrix<double>::eigRight(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi, 
    Matrix<double>& vec)
  {
// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(numRows()*numColumns());
    unsigned count = 0;
    for (int j=getLower(1); j<getUpper(1); ++j)
      for (int i=getLower(0); i<getUpper(0); ++i)
        A[count++] = this->operator()(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = 1;
    LDVR = vec.numRows();

    LDA = numRows();
    N = LDA;
    JOBVR = 'V';
    JOBVL = 'N';
    LWORK = 5*N;
    std::vector<double> WORK(LWORK);
     
// call LAPACK routine to compute eigenvalues and eigenvectors
    dgeev_(&JOBVL, &JOBVR, &N, &A[0], &LDA,
      &evr[evr.getLower(0)],
      &evi[evi.getLower(0)],
      0, &LDVL,
      &vec(vec.getLower(0), vec.getLower(1)), &LDVR,
      &WORK[0], &LWORK,
      &INFO);

// check if eigensolver worked
    if (INFO!=0)
      throw Lucee::Except("Matrix::eigRight: Eigenvalue solver failed");
  }

  template <>
  void 
  Matrix<double>::eigLeft(Lucee::Vector<double>& evr, Lucee::Vector<double>& evi,
    Matrix<double>& vec)
  {
// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(numRows()*numColumns());
    unsigned count = 0;
    for (int j=getLower(1); j<getUpper(1); ++j)
      for (int i=getLower(0); i<getUpper(0); ++i)
        A[count++] = this->operator()(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = vec.numRows();
    LDVR = 1;

    LDA = numRows();
    N = LDA;
    JOBVR = 'N';
    JOBVL = 'V';
    LWORK = 5*N;
    std::vector<double> WORK(LWORK);
     
// call LAPACK routine to compute eigenvalues and eigenvectors
    dgeev_(&JOBVL, &JOBVR, &N, &A[0], &LDA,
      &evr[evr.getLower(0)],
      &evi[evi.getLower(0)],
      &vec(vec.getLower(0), vec.getLower(1)), &LDVL,
      0, &LDVR,
      &WORK[0], &LWORK,
      &INFO);

// check if eigensolver worked
    if (INFO!=0)
      throw Lucee::Except("Matrix::eigLeft: Eigenvalue solver failed");
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
