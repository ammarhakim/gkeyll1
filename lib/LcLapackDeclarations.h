/**
 * @file	LcLapackDeclarations.h
 *
 * @brief	Declarations for pulling in LAPACK symbols into Lucee.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_LAPACK_DECLARATIONS_H
#define LC_LAPACK_DECLARATIONS_H

extern "C"
{
/** DGEEV */
// *  Purpose
// *  =======
// *
// *  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
// *  eigenvalues and, optionally, the left and/or right eigenvectors.
// *
// *  The right eigenvector v(j) of A satisfies
// *                   A * v(j) = lambda(j) * v(j)
// *  where lambda(j) is its eigenvalue.
// *  The left eigenvector u(j) of A satisfies
// *                u(j)**H * A = lambda(j) * u(j)**H
// *  where u(j)**H denotes the conjugate transpose of u(j).
// *
// *  The computed eigenvectors are normalized to have Euclidean norm
// *  equal to 1 and largest component real.
// *
// *  Arguments
// *  =========
// *
// *  JOBVL   (input) CHARACTER*1
// *          = 'N': left eigenvectors of A are not computed;
// *          = 'V': left eigenvectors of A are computed.
// *
// *  JOBVR   (input) CHARACTER*1
// *          = 'N': right eigenvectors of A are not computed;
// *          = 'V': right eigenvectors of A are computed.
// *
// *  N       (input) INTEGER
// *          The order of the matrix A. N >= 0.
// *
// *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
// *          On entry, the N-by-N matrix A.
// *          On exit, A has been overwritten.
// *
// *  LDA     (input) INTEGER
// *          The leading dimension of the array A.  LDA >= max(1,N).
// *
// *  WR      (output) DOUBLE PRECISION array, dimension (N)
// *  WI      (output) DOUBLE PRECISION array, dimension (N)
// *          WR and WI contain the real and imaginary parts,
// *          respectively, of the computed eigenvalues.  Complex
// *          conjugate pairs of eigenvalues appear consecutively
// *          with the eigenvalue having the positive imaginary part
// *          first.
// *
// *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
// *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
// *          after another in the columns of VL, in the same order
// *          as their eigenvalues.
// *          If JOBVL = 'N', VL is not referenced.
// *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
// *          the j-th column of VL.
// *          If the j-th and (j+1)-st eigenvalues form a complex
// *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
// *          u(j+1) = VL(:,j) - i*VL(:,j+1).
// *
// *  LDVL    (input) INTEGER
// *          The leading dimension of the array VL.  LDVL >= 1; if
// *          JOBVL = 'V', LDVL >= N.
// *
// *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
// *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
// *          after another in the columns of VR, in the same order
// *          as their eigenvalues.
// *          If JOBVR = 'N', VR is not referenced.
// *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
// *          the j-th column of VR.
// *          If the j-th and (j+1)-st eigenvalues form a complex
// *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
// *          v(j+1) = VR(:,j) - i*VR(:,j+1).
// *
// *  LDVR    (input) INTEGER
// *          The leading dimension of the array VR.  LDVR >= 1; if
// *          JOBVR = 'V', LDVR >= N.
// *
// *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
// *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
// *
// *  LWORK   (input) INTEGER
// *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
// *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
// *          performance, LWORK must generally be larger.
// *
// *          If LWORK = -1, then a workspace query is assumed; the routine
// *          only calculates the optimal size of the WORK array, returns
// *          this value as the first entry of the WORK array, and no error
// *          message related to LWORK is issued by XERBLA.
// *
// *  INFO    (output) INTEGER
// *          = 0:  successful exit
// *          < 0:  if INFO = -i, the i-th argument had an illegal value.
// *          > 0:  if INFO = i, the QR algorithm failed to compute all the
// *                eigenvalues, and no eigenvectors have been computed;
// *                elements i+1:N of WR and WI contain eigenvalues which
// *                have converged.
// *
// *  =====================================================================
    void dgeev_(char *jobvl, char *jobvr,
      int *n, double *a, int *lda,
      double *wr, double *wi,
      double *vl, int *ldvl,
      double *vr, int *ldvr,
      double *work, int *lwork, int *info);

/** DGEMM */
// *
// *  Purpose
// *  =======
// *
// *  DGEMM  performs one of the matrix-matrix operations
// *
// *     C := alpha*op( A )*op( B ) + beta*C,
// *
// *  where  op( X ) is one of
// *
// *     op( X ) = X   or   op( X ) = X',
// *
// *  alpha and beta are scalars, and A, B and C are matrices, with op( A )
// *  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
// *
// *  Arguments
// *  ==========
// *
// *  TRANSA - CHARACTER*1.
// *           On entry, TRANSA specifies the form of op( A ) to be used in
// *           the matrix multiplication as follows:
// *
// *              TRANSA = 'N' or 'n',  op( A ) = A.
// *
// *              TRANSA = 'T' or 't',  op( A ) = A'.
// *
// *              TRANSA = 'C' or 'c',  op( A ) = A'.
// *
// *           Unchanged on exit.
// *
// *  TRANSB - CHARACTER*1.
// *           On entry, TRANSB specifies the form of op( B ) to be used in
// *           the matrix multiplication as follows:
// *
// *              TRANSB = 'N' or 'n',  op( B ) = B.
// *
// *              TRANSB = 'T' or 't',  op( B ) = B'.
// *
// *              TRANSB = 'C' or 'c',  op( B ) = B'.
// *
// *           Unchanged on exit.
// *
// *  M      - INTEGER.
// *           On entry,  M  specifies  the number  of rows  of the  matrix
// *           op( A )  and of the  matrix  C.  M  must  be at least  zero.
// *           Unchanged on exit.
// *
// *  N      - INTEGER.
// *           On entry,  N  specifies the number  of columns of the matrix
// *           op( B ) and the number of columns of the matrix C. N must be
// *           at least zero.
// *           Unchanged on exit.
// *
// *  K      - INTEGER.
// *           On entry,  K  specifies  the number of columns of the matrix
// *           op( A ) and the number of rows of the matrix op( B ). K must
// *           be at least  zero.
// *           Unchanged on exit.
// *
// *  ALPHA  - DOUBLE PRECISION.
// *           On entry, ALPHA specifies the scalar alpha.
// *           Unchanged on exit.
// *
// *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
// *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
// *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
// *           part of the array  A  must contain the matrix  A,  otherwise
// *           the leading  k by m  part of the array  A  must contain  the
// *           matrix A.
// *           Unchanged on exit.
// *
// *  LDA    - INTEGER.
// *           On entry, LDA specifies the first dimension of A as declared
// *           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
// *           LDA must be at least  max( 1, m ), otherwise  LDA must be at
// *           least  max( 1, k ).
// *           Unchanged on exit.
// *
// *  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
// *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
// *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
// *           part of the array  B  must contain the matrix  B,  otherwise
// *           the leading  n by k  part of the array  B  must contain  the
// *           matrix B.
// *           Unchanged on exit.
// *
// *  LDB    - INTEGER.
// *           On entry, LDB specifies the first dimension of B as declared
// *           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
// *           LDB must be at least  max( 1, k ), otherwise  LDB must be at
// *           least  max( 1, n ).
// *           Unchanged on exit.
// *
// *  BETA   - DOUBLE PRECISION.
// *           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
// *           supplied as zero then C need not be set on input.
// *           Unchanged on exit.
// *
// *  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
// *           Before entry, the leading  m by n  part of the array  C must
// *           contain the matrix  C,  except when  beta  is zero, in which
// *           case C need not be set on entry.
// *           On exit, the array  C  is overwritten by the  m by n  matrix
// *           ( alpha*op( A )*op( B ) + beta*C ).
// *
// *  LDC    - INTEGER.
// *           On entry, LDC specifies the first dimension of C as declared
// *           in  the  calling  (sub)  program.   LDC  must  be  at  least
// *           max( 1, m ).
// *           Unchanged on exit.
// *
    void dgemm_(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
      double *ALPHA, double A[], int *LDA,
      double B[], int *LDB, double *beta,
      double C[], int *LDC);

// *
// *  Purpose
// *  =======
// *
// *  DGEMV  performs one of the matrix-vector operations
// *
// *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
// *
// *  where alpha and beta are scalars, x and y are vectors and A is an
// *  m by n matrix.
// *
// *  Arguments
// *  ==========
// *
// *  TRANS  - CHARACTER*1.
// *           On entry, TRANS specifies the operation to be performed as
// *           follows:
// *
// *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
// *
// *              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
// *
// *              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
// *
// *           Unchanged on exit.
// *
// *  M      - INTEGER.
// *           On entry, M specifies the number of rows of the matrix A.
// *           M must be at least zero.
// *           Unchanged on exit.
// *
// *  N      - INTEGER.
// *           On entry, N specifies the number of columns of the matrix A.
// *           N must be at least zero.
// *           Unchanged on exit.
// *
// *  ALPHA  - DOUBLE PRECISION.
// *           On entry, ALPHA specifies the scalar alpha.
// *           Unchanged on exit.
// *
// *  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
// *           Before entry, the leading m by n part of the array A must
// *           contain the matrix of coefficients.
// *           Unchanged on exit.
// *
// *  LDA    - INTEGER.
// *           On entry, LDA specifies the first dimension of A as declared
// *           in the calling (sub) program. LDA must be at least
// *           max( 1, m ).
// *           Unchanged on exit.
// *
// *  X      - DOUBLE PRECISION array of DIMENSION at least
// *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
// *           and at least
// *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
// *           Before entry, the incremented array X must contain the
// *           vector x.
// *           Unchanged on exit.
// *
// *  INCX   - INTEGER.
// *           On entry, INCX specifies the increment for the elements of
// *           X. INCX must not be zero.
// *           Unchanged on exit.
// *
// *  BETA   - DOUBLE PRECISION.
// *           On entry, BETA specifies the scalar beta. When BETA is
// *           supplied as zero then Y need not be set on input.
// *           Unchanged on exit.
// *
// *  Y      - DOUBLE PRECISION array of DIMENSION at least
// *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
// *           and at least
// *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
// *           Before entry with BETA non-zero, the incremented array Y
// *           must contain the vector y. On exit, Y is overwritten by the
// *           updated vector y.
// *
// *  INCY   - INTEGER.
// *           On entry, INCY specifies the increment for the elements of
// *           Y. INCY must not be zero.
// *           Unchanged on exit.
    
    void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double A[], int *LDA, 
      double X[], int *INCX, double *BETA, double Y[], int *INCY);
}

#endif // LC_LAPACK_DECLARATIONS_H
