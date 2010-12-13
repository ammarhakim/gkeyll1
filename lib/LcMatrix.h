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
#include <LcColMajorIndexer.h>
#include <LcVector.h>

namespace Lucee
{
/**
 * Class to represent matrices.
 */
  template <typename T>
  class Matrix : public Lucee::Array<2, T, Lucee::ColMajorIndexer >
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
 * Creat a new matrix reusing supplied memory. This constructor can
 * not check if the data provided is sufficient to create the matrix
 * and hence it is the user's reposibility to ensure this.
 *
 * @param row Rows in matrix.
 * @param col Columns in matrix.
 * @param data Memory to reuse. This should be atleast row X col X sizeof(T).
 */
      Matrix(unsigned row, unsigned col, T *data);

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
 * Duplicate this matrix.
 *
 * @return Copy of this matrix.
 */
      Matrix<T> duplicate() const;

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
 * Return a row of the martrix.
 *
 * @param row Row to return.
 * @return Row of the matrix.
 */
      Lucee::Vector<T> getRow(unsigned row);

/**
 * Return a column of the martrix.
 *
 * @param col Column to return.
 * @return Column of the matrix.
 */
      Lucee::Vector<T> getCol(unsigned col);

/**
 * Scale all rows of the matrix with by multipling with elements of
 * supplied vector. I.e. M[i,j] <- M[i,j]*fv[i];
 *
 * @param fv Vector with same number of elements as number of rows.
 */
      void scaleRows(const Lucee::Vector<double>& fv);

/**
 * Scale all columns of the matrix with by multipling with elements of
 * supplied vector. I.e. M[i,j] <- M[i,j]*fv[j];
 *
 * @param fv Vector with same number of elements as number of columns.
 */
      void scaleCols(const Lucee::Vector<double>& fv);

/**
 * Get view into matrix. The returned matrix shares data with this
 * matrix.
 * 
 * @param lr Lower row index into view.
 * @param ur Upper row index into view.
 * @param lc Lower column index into view.
 * @param uc Upper column index index view.
 * @return view matrix.
 */
      Lucee::Matrix<T> getView(int lr, int ur, int lc, int uc);

/**
 * Is this matrix a transpose of another one?
 *
 * @return true if this matrix is transpose of another one.
 */
      bool isTranspose() const;

/**
 * Is this matrix square?
 *
 * @return true if this matrix is square, false otherwise.
 */
      bool isSquare() const { return numColumns() == numRows(); }

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

/**
 * Create matrix from given array.
 *
 * @param arr Array to create from.
 */
      Matrix(Lucee::Array<2, T, Lucee::ColMajorIndexer>& arr);
  };
}

#endif //  LC_MATRIX_H
