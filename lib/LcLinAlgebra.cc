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
  Lucee::Matrix<double>&
  accumulate(double beta, Lucee::Matrix<double>& C,
    double alpha, const Lucee::Matrix<double>& A, const Lucee::Matrix<double>& B)
  {
    unsigned crows = C.numRows();
    unsigned ccols = C.numColumns();
    unsigned arows = A.numRows();
    unsigned acols = A.numColumns();
    unsigned brows = B.numRows();
    unsigned bcols = B.numColumns();

    char TRANSA, TRANSB;
    TRANSA = 'N'; TRANSB = 'N'; // by default do not transpose A and B

// Check that shapes of matrices are consistent. One needs to be
// careful as matrices may be transposed
    if (A.isTranspose())
    {
      arows = A.numColumns();
      acols = A.numRows();
      TRANSA = 'T';
    }
    if (B.isTranspose())
    {
      brows = B.numColumns();
      bcols = B.numRows();
      TRANSB = 'T';
    }

    if ((arows != crows) || (bcols != ccols) || (acols != brows))
      throw Lucee::Except("Lucee::accumulate: Inconsistent shape of matrices");

// Check if input and output matrices are contiguous.
    Matrix<double> Cdup(C);
    if (C.isContiguous() == false)
      Cdup = C.duplicate(); // not, so allocate fresh matrix
    Matrix<double> Adup(A);
    if (A.isContiguous() == false)
      Adup = A.duplicate(); // not, so allocate fresh matrix
    Matrix<double> Bdup(B);
    if (B.isContiguous() == false)
      Bdup = B.duplicate(); // not, so allocate fresh matrix

    return C;
  }


  void 
  eig(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi)
  {
    if (mat.isSquare()==false)
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
      &evr.first(),
      &evi.first(),
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
    if (mat.isSquare()==false)
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
      &evr.first(),
      &evi.first(),
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
    if (mat.isSquare()==false)
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
      &evr.first(),
      &evi.first(),
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
    if (mat.isSquare()==false)
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
      &evr.first(),
      &evi.first(),
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
