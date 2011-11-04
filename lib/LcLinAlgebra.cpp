/**
 * @file	LcLinAlgebra.cpp
 *
 * @brief	Linear algebra functions. Mainly wrappers around LAPACK.
 */

// lucee includes
#include <LcExcept.h>
#include <LcLinAlgebra.h>
#include <LcLapackDeclarations.h>

// std includes
#include <vector>
#include <cstdlib>

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

// check that shapes of matrices are consistent
    if ((arows != crows) || (bcols != ccols) || (acols != brows))
      throw Lucee::Except("Lucee::accumulate: Inconsistent shape of matrices");

// stuff needed by BLAS routine
    int M, N, K, LDA, LDB, LDC;
    char TRANSA, TRANSB;
    TRANSA = 'N'; TRANSB = 'N'; // by default do not transpose A and B

    M = arows;
    N = bcols;
    K = acols;
    LDA = arows;
    LDB = brows;
    LDC = crows;

// determine job flags and sizes to pass to BLAS
    if (A.isTranspose())
    {
      TRANSA = 'T';
      LDA = acols; // underlying data is not really transposed
    }
    if (B.isTranspose())
    {
      TRANSB = 'T';
      LDB = bcols; // underlying data is not really transposed
    }

// check if input and output matrices are contiguous
    Matrix<double> Cdup(C);
    if (C.isContiguous() == false)
      Cdup = C.duplicate(); // not, so allocate fresh matrix
    Matrix<double> Adup(A);
    if (A.isContiguous() == false)
      Adup = A.duplicate(); // not, so allocate fresh matrix
    Matrix<double> Bdup(B);
    if (B.isContiguous() == false)
      Bdup = B.duplicate(); // not, so allocate fresh matrix

// call BLAS routine to do the multiplication
    dgemm_(&TRANSA, &TRANSB, &M, &N, &K,
      &alpha, &Adup.first(), &LDA, &Bdup.first(), &LDB,
      &beta, &Cdup.first(), &LDC);

    if (C.isContiguous() == false)
    {
// C was not contiguous, so copy data from Cdup to C
      for (int i=C.getLower(0); i<C.getUpper(0); ++i)
        for (int j=C.getLower(1); j<C.getUpper(1); ++j)
          C(i,j) = Cdup(i,j);
    }
    return C;
  }

  Lucee::Matrix<double>& accumulate(Lucee::Matrix<double>& C,
    const Lucee::Matrix<double>& A, const Lucee::Matrix<double>& B)
  {
    return accumulate(0.0, C, 1.0, A, B);
  }

  Lucee::Vector<double>&
  accumulate(double beta, Lucee::Vector<double>& y,
    double alpha, const Lucee::Matrix<double>& A, const Lucee::Vector<double>& x)
  {
    unsigned ylen = y.getLength();
    unsigned arows = A.numRows();
    unsigned acols = A.numColumns();
    unsigned xlen = x.getLength();

// check that shapes of matrices and vectors are consistent
    if ((acols != xlen) || (ylen != arows))
      throw Lucee::Except("Lucee::accumulate: Inconsistent shape of matrix and vector.");

    char TRANS;
    int M, N, LDA, INCX, INCY;
    TRANS = 'N';

    M = arows;
    N = acols;
    LDA = arows;
    INCX = 1;
    INCY = 1;

// determine flag to pass to BLAS
    if (A.isTranspose())
    {
      TRANS = 'T';
      LDA = acols;
    }

// check if input and output matrices/vectors are contigous
    Lucee::Matrix<double> Adup(A);
    if (A.isContiguous() == false)
      Adup = A.duplicate(); // not, so allocate fresh matrix
    Lucee::Vector<double> xdup(x);
    if (x.isContiguous() == false)
      xdup = x.duplicate(); // not, so allocate fresh vector
    Lucee::Vector<double> ydup(y);
    if (y.isContiguous() == false)
      ydup = y.duplicate(); // not, so allocate fresh vector

// call BLAS routine to do the multiplication    
    dgemv_(&TRANS, &M, &N, &alpha, &Adup.first(), &LDA,
      &xdup.first(), &INCX, &beta, &ydup.first(), &INCY);

    if (y.isContiguous() == false)
    {
// y was not contigous so copy data from ydup to y
      for (int i=y.getLower(0); i<y.getUpper(0); ++i)
        y[i] = ydup[i];
    }
    return y;
  }

  Lucee::Vector<double>& accumulate(Lucee::Vector<double>& y,
    const Lucee::Matrix<double>& A, const Lucee::Vector<double>& x)
  {
    return accumulate(0.0, y, 1.0, A, x);
  }

  Lucee::Matrix<double>&
  accumulate(Lucee::Matrix<double>& A, double alpha, 
    const Lucee::Vector<double>& x, const Lucee::Vector<double>& y)
  {
    unsigned arows = A.numRows();
    unsigned acols = A.numColumns();
    unsigned xlen = x.getLength();
    unsigned ylen = y.getLength();

// check that shapes of matrices and vectors are consistent
    if ((arows != xlen) || (acols != ylen))
      throw Lucee::Except("Lucee::accumulate: Inconsistent shape of matrices");

    int M, N, INCX, INCY, LDA;
    M = arows;
    N = acols;
    INCX = 1;
    INCY = 1;
    LDA = arows;

    if (A.isTranspose())
      LDA = acols;

// check if input and output matrices/vectors are contigous
    Lucee::Matrix<double> Adup(A);
    if (A.isContiguous() == false)
      Adup = A.duplicate(); // not, so allocate fresh matrix
    Lucee::Vector<double> xdup(x);
    if (x.isContiguous() == false)
      xdup = x.duplicate(); // not, so allocate fresh vector
    Lucee::Vector<double> ydup(y);
    if (y.isContiguous() == false)
      ydup = y.duplicate(); // not, so allocate fresh vector

// call BLAS routine to do the multiplication
    dger_(&M, &N, &alpha, &xdup.first(), &INCX, &ydup.first(), &INCY, &Adup.first(), &LDA);

    if (A.isContiguous() == false)
    {
// A is not contiguous, so copy data from Adup to A
      for (int i=A.getLower(0); i<A.getUpper(0); ++i)
        for (int j=A.getLower(1); j<A.getUpper(1); ++j)
          A(i,j) = Adup(i,j);
    }
    return A;
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
  solve(Lucee::Matrix<double>& A, Lucee::Matrix<double>& B)
  {
// check pre-conditions
    if (A.isSquare() == false)
      throw Lucee::Except("Lucee::solve: LHS matrix must be square.");
    if (A.isContiguous() == false)
      throw Lucee::Except("Lucee::solve: LHS matrix must be contiguous.");
    if (B.isContiguous() == false)
      throw Lucee::Except("Lucee::solve: RHS matrix must be contiguous.");

    int INFO, LDA, LDB, M, N, NRHS;
    int *IPIV;
    char TRANS;

    M = N = LDA = LDB = A.getShape(0);
// allocate memory for permutation array
    IPIV = (int*) malloc( sizeof(int)*M );

// compute LU factorization
    dgetrf_(&M, &N, &A.first(), &LDA, IPIV, &INFO);
    if(INFO != 0)
      throw Lucee::Except("Lucee::solve: Unable to compute LU factorization");

// solve systems
    TRANS = 'N';
    NRHS = B.getShape(1);
// make call
    dgetrs_(&TRANS, &N, &NRHS, &A.first(), &LDA, IPIV, &B.first(), &LDB, &INFO);
// check if solution obtained okay
    if (INFO != 0)
      throw Lucee::Except("Lucee::solve: Failure in LAPACK linear solver");

    free(IPIV);
  }
}
