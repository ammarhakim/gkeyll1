/**
 * @file	LcLinAlgebra.h
 *
 * @brief	Linear algebra functions. Mainly wrappers around LAPACK.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_LIN_ALGEBRA_H
#define LC_LIN_ALGEBRA_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcMatrix.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Computes matrix-matrix product C = alpha*A*B + beta*C.
 *
 * @param beta Coefficient in front of C.
 * @param C output/input matrix.
 * @param alpha Coefficient in from on A*B.
 * @param A Input matrix appearing in A*B.
 * @param B Input matrix appearing in A*B.
 * @return Reference to updated matrix C.
 */
  Lucee::Matrix<double>& accumulate(double beta, Lucee::Matrix<double>& C,
    double alpha, const Lucee::Matrix<double>& A, const Lucee::Matrix<double>& B);

/**
 * Computes matrix-vector product y = alpha*A*x + beta*y, where A is a
 * matrix and x and y are vectors.
 *
 * @param beta Coefficient in front of y.
 * @param y output/input vector.
 * @param alpha Coefficient in from on A*x.
 * @param A Input matrix appearing in A*x.
 * @param x Input matrix appearing in A*x.
 * @return Reference to updated vector y.
 */
  Lucee::Vector<double>& accumulate(double beta, Lucee::Vector<double>& y,
    double alpha, const Lucee::Matrix<double>& A, const Lucee::Vector<double>& x);

/**
 * Computes vector-vector outer product A = alpha*x*y' + A, where A is
 * a matrix and x and y are vectors.
 *
 * @param A input/output matrix
 * @param alpha Coefficient in from on x*y'.
 * @param x Input matrix appearing in x*y'.
 * @param y Input matrix appearing in x*y'.
 * @return Reference to updated matrix A.
 */
  Lucee::Matrix<double>& accumulate(Lucee::Matrix<double>& A,
    double alpha, const Lucee::Vector<double>& x, const Lucee::Vector<double>& y);

/**
 * Computes eigenvalues of a matrix. Matrix must be square or an
 * exception will be thrown.
 *
 * @param mat Matrix. Must be square.
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 */
  void eig(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi);

/**
 * Computes eigenvalues and left, right eigenvectors of matrix. Matrix
 * must be square or an exception will be thrown.
 *
 * @param mat Matrix. Must be square.
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vecl Left eigenvectors of matrix.
 * @param vecr Right eigenvectors of matrix.
 */
  void eig(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi, Matrix<double>& vecl, Matrix<double>& vecr);

/**
 * Computes eigenvalues and right eigenvectors of matrix. Matrix must
 * be square or an exception will be thrown.
 *
 * @param mat Matrix. Must be square.
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vec Left eigenvectors of matrix.
 */
  void eigRight(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi, Lucee::Matrix<double>& vec);

/**
 * Computes eigenvalues and left eigenvectors of matrix. Matrix must
 * be square or an exception will be thrown.
 *
 * @param mat Matrix. Must be square.
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vec Left eigenvectors of matrix.
 */
  void eigLeft(const Lucee::Matrix<double>& mat, Lucee::Vector<double>& evr, 
    Lucee::Vector<double>& evi, Lucee::Matrix<double>& vec);

/**
 * Solve system of linear equations Ax = b, where A is 'mat' and 'b'
 * are the column vectors of the 'rhs' matrix. On output 'rhs'
 * contains the corresponding solution vectors.
 *
 * @param mat Matrix in linear system Ax = b.  
 * @param rhs Column vectors of 'rhs' are 'b'. On output 'rhs' is
 *   replaced by the corresponding solution vectors.
 */
  void solve(const Lucee::Matrix<double>& mat, Lucee::Matrix<double>& rhs);
}

#endif // LC_LIN_ALGEBRA_H
