/**
 * @file	LcLinAlgebra.cc
 *
 * @brief	Linear algebra functions. Mainly wrappers around LAPACK.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

// lucee includes
#include <LcExcept.h>
#include <LcLinAlgebra.h>
#include <LcLapackDeclarations.h>

// std includes
#include <vector>

namespace Lucee
{
  template <typename T>
  static bool checkIfMatrixIsSquare(const Lucee::Matrix<T>& mat)
  {
    return mat.numRows() == mat.numColumns();
  }

  void 
  eig(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi)
  {
    if (checkIfMatrixIsSquare(mat)==false)
      throw Lucee::Except("Lucee:eig: Matrix must be square.");

// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(mat.numRows()*mat.numColumns());
    unsigned count = 0;
    for (int j=mat.getLower(1); j<mat.getUpper(1); ++j)
      for (int i=mat.getLower(0); i<mat.getUpper(0); ++i)
        A[count++] = mat(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = 1;
    LDVR = 1;

    LDA = mat.numRows();
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
      throw Lucee::Except("Lucee::eig: Eigenvalue solver failed");
  }

  void
  eig(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi, Matrix<double>& vecl, Matrix<double>& vecr)
  {
    if (checkIfMatrixIsSquare(mat)==false)
      throw Lucee::Except("Lucee:eig: Matrix must be square.");

// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(mat.numRows()*mat.numColumns());
    unsigned count = 0;
    for (int j=mat.getLower(1); j<mat.getUpper(1); ++j)
      for (int i=mat.getLower(0); i<mat.getUpper(0); ++i)
        A[count++] = mat(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = vecl.numRows();
    LDVR = vecl.numRows();

    LDA = mat.numRows();
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

  void
  eigRight(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi, Matrix<double>& vec)
  {
    if (checkIfMatrixIsSquare(mat)==false)
      throw Lucee::Except("Lucee:eig: Matrix must be square.");

// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(mat.numRows()*mat.numColumns());
    unsigned count = 0;
    for (int j=mat.getLower(1); j<mat.getUpper(1); ++j)
      for (int i=mat.getLower(0); i<mat.getUpper(0); ++i)
        A[count++] = mat(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = 1;
    LDVR = vec.numRows();

    LDA = mat.numRows();
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

  void
  eigLeft(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr,
    Lucee::Vector<double>& evi, Matrix<double>& vec)
  {
    if (checkIfMatrixIsSquare(mat)==false)
      throw Lucee::Except("Lucee:eig: Matrix must be square.");

// copy data from matrix into temporary array. Use column major order.
    std::vector<double> A(mat.numRows()*mat.numColumns());
    unsigned count = 0;
    for (int j=mat.getLower(1); j<mat.getUpper(1); ++j)
      for (int i=mat.getLower(0); i<mat.getUpper(0); ++i)
        A[count++] = mat(i,j);

    char JOBVR, JOBVL;
    int INFO, LWORK, LDA, N, LDVL, LDVR;
    LDVL = vec.numRows();
    LDVR = 1;

    LDA = mat.numRows();
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

  void
  solve(const Lucee::Matrix<double>& mat, Lucee::Matrix<double>& rhs)
  {
  }
}
