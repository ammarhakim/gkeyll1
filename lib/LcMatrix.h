/**
 * @file	LcMatrix.h
 *
 * @brief	Matrix class.
 *
 * @version	$Id$
 *
 * Copyright &copy; 2008-2009, Ammar Hakim.
 */

#ifndef LC_MATRIX_H
#define LC_MATRIX_H

// config stuff
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

// lucee includes
#include <LcArray.h>

namespace Lucee
{
  template <typename T>
  class Matrix : public Lucee::Array<2, T, Lucee::ColMajorIndexer<2> >
  {
    public:
/**
 * Construct matrix with specified number of rows and columns.
 *
 * @param row Rows in matrix.
 * @param col Columns in matrix.
 */      
      Matrix(unsigned row, unsigned col);

/**
 * Construct matrix with specified shape.
 *
 * @param shape Shape of the matrix.
 */      
      Matrix(unsigned shape[2]);

/**
 * Construct matrix with specified shape and start indices.
 *
 * @param shape Shape of the matrix.
 * @param start Start indices for the matrix.
 */
      Matrix(unsigned shape[2], int start[2]);

/**
 * Returns number of rows in matrix.
 *
 * @return number of rows.
 */
      unsigned numRows() const { return this->template getShape(0); }

/**
 * Returns number of columns in matrix.
 *
 * @return number of columns.
 */
      unsigned numColumns() const { return this->template getShape(1); }

/**
 * Computes eigenvalues of the matrix. Matrix must be square or an
 * exception will be thrown.
 *
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 */
      void eig(Lucee::Array<1, T>& evr, Lucee::Array<1, T>& evi);

/**
 * Solves the linear system of equations Ax = b, where A is this
 * matrix and b are column vectors of the 'rhs' matrix. On output the
 * 'rhs' matrix is replaced by the solution vectors. An LU
 * decomposition algorithm is used.
 *
 * @param rhs On input columns contain RHS and on output solution.
 */
      void solve(Lucee::Matrix<T>& rhs);
  };
}

#endif //  LC_MATRIX_H
