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
#include <LcVector.h>

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
 * Creat a new matrix from input matrix.
 *
 * @param mat Matrix to copy from
 */
      Matrix(const Matrix<T>& mat);

/**
 * Copy input matrix.
 *
 * @param mat Matrix to copy.
 * @return Reference to this matrix.
 */
      Matrix<T>& operator=(const Matrix<T>& mat);

/**
 * Assign all elements in matrix to specified value.
 *
 * @param val Value to assign.
 * @return Reference to this matrix.
 */
      Matrix<T>& operator=(const T& val);

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
 * Return the transpose the matrix. No data is actually allocated and
 * the transpose matrix and shares data with this matrix.
 *
 * @return transpose of this matrix.
 */
      Matrix<T> transpose() const;

/**
 * Computes eigenvalues of the matrix. Matrix must be square or an
 * exception will be thrown.
 *
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 */
      void eig(Lucee::Vector<T>& evr, Lucee::Vector<T>& evi);

/**
 * Computes eigenvalues and left, right eigenvectors of matrix. Matrix
 * must be square or an exception will be thrown.
 *
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vecl Left eigenvectors of matrix.
 * @param vecr Right eigenvectors of matrix.
 */
      void eig(Lucee::Vector<T>& evr, Lucee::Vector<T>& evi, 
        Matrix<T>& vecl, Matrix<T>& vecr);

/**
 * Computes eigenvalues and right eigenvectors of matrix. Matrix must
 * be square or an exception will be thrown.
 *
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vec Left eigenvectors of matrix.
 */
      void eigRight(Lucee::Vector<T>& evr, Lucee::Vector<T>& evi, 
        Matrix<T>& vec);

/**
 * Computes eigenvalues and left eigenvectors of matrix. Matrix must
 * be square or an exception will be thrown.
 *
 * @param evr Real part of eigenvalues.
 * @param evi Imaginary part of eigenvalues.
 * @param vec Left eigenvectors of matrix.
 */
      void eigLeft(Lucee::Vector<T>& evr, Lucee::Vector<T>& evi,
        Matrix<T>& vec);

/**
 * Solves the linear system of equations Ax = b, where A is this
 * matrix and b are column vectors of the 'rhs' matrix. On output the
 * 'rhs' matrix is replaced by the solution vectors. An LU
 * decomposition algorithm is used.
 *
 * @param rhs On input columns contain RHS and on output solution.
 */
      void solve(Lucee::Matrix<T>& rhs);

    private:
/** Matrix traits stored as bit values */
      unsigned traits;
  };

}

#endif //  LC_MATRIX_H
