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
