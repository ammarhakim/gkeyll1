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

    private:
/** Matrix traits stored as bit values */
      unsigned traits;
  };
}

#endif //  LC_MATRIX_H
